% =========================================================================
% ComponentWiseHeldOutCorrelation.m
%
% Purpose:
%   Evaluate the component-wise held-out performance of CCA, G-FWCCA,
%   LG-FWCCA, and GL-FWCCA for the resting-state fMRI analysis.
%
%   For each subject and each default mode network (DMN) ROI, this script
%   computes the absolute Pearson correlation between:
%
%       (1) the held-out seed-correlation map, and
%       (2) each held-out canonical spatial map.
%
%   The analysis is performed separately for the first four canonical
%   components. Group-level results are summarized as the mean +/- standard
%   error (SE) across the five subjects.
%
% ROIs:
%   AnG  : Angular Gyrus
%   PCgG : Posterior Cingulate Gyrus
%   PCu  : Precuneus
%   SFG  : Superior Frontal Gyrus
%
% Methods:
%   CCA
%   G-FWCCA
%   LG-FWCCA
%   GL-FWCCA
%
% Input:
%   ../03_ModelTesting/Results/rsfMRI_TestingResults.mat
%
% Required variable:
%   AllTestMaps
%
% Output:
%   Results/ComponentWiseCorrResults.mat
%
%   Results/ComponentWiseCorrelationFigures/
%       AnG_ComponentWise_HeldOutSeedCorrelation.pdf/.png
%       PCgG_ComponentWise_HeldOutSeedCorrelation.pdf/.png
%       PCu_ComponentWise_HeldOutSeedCorrelation.pdf/.png
%       SFG_ComponentWise_HeldOutSeedCorrelation.pdf/.png
%
% Note:
%   This script performs only the component-wise held-out correlation
%   analysis. The representative-subject spatial-overlap analysis is
%   performed separately.
%
% =========================================================================

clear;
clc;
close all;


%% ========================================================================
% Project paths
% =========================================================================

% Directory containing this script:
%   openNeuro_ds005747/04_ModelEvaluation
scriptDir = fileparts(mfilename('fullpath'));

% Dataset project directory:
%   openNeuro_ds005747
datasetRoot = fileparts(scriptDir);

% Step 03 testing-results directory
testingResultsDir = fullfile( ...
    datasetRoot, ...
    '03_ModelTesting', ...
    'Results');

% Input generated in Step 03
resultFile = fullfile( ...
    testingResultsDir, ...
    'rsfMRI_TestingResults.mat');

% Output directory for Step 04
outputDir = fullfile(scriptDir, 'Results');

% Figure output directory
figureDir = fullfile( ...
    outputDir, ...
    'ComponentWiseCorrelationFigures');


%% ========================================================================
% Check required files/directories
% =========================================================================

if ~isfile(resultFile)
    error(['Step 03 testing-result file not found:\n%s\n\n' ...
           'Run rsfMRI_Testing.m first.'], ...
           resultFile);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end


%% ========================================================================
% Load held-out testing maps
% =========================================================================

fprintf('============================================================\n');
fprintf('Component-wise Held-Out Correlation Analysis\n');
fprintf('============================================================\n\n');

fprintf('Loading Step 03 testing results:\n%s\n\n', resultFile);

S = load(resultFile, 'AllTestMaps');

if ~isfield(S, 'AllTestMaps')
    error('Variable "AllTestMaps" was not found in:\n%s', resultFile);
end

AllTestMaps = S.AllTestMaps;


%% ========================================================================
% Analysis settings
% =========================================================================

roi_list = { ...
    'AnG', ...
    'PCgG', ...
    'PCu', ...
    'SFG'};

methodNames = { ...
    'CCA', ...
    'G-FWCCA', ...
    'LG-FWCCA', ...
    'GL-FWCCA'};

fieldNames = { ...
    'CCA', ...
    'G_FWCCA', ...
    'LG_FWCCA', ...
    'GL_FWCCA'};

markerList = { ...
    'o', ...     % CCA
    '^', ...     % G-FWCCA
    's', ...     % LG-FWCCA
    'd'};        % GL-FWCCA

Kcomp = 4;

xBase = 1:Kcomp;

% Horizontal shifts used only for visualization
xShift = [-0.22, -0.05, 0.10, 0.22];


%% ========================================================================
% Initialise output structure
% =========================================================================

ComponentWiseCorrResults = struct();


%% ========================================================================
% Loop over ROIs
% =========================================================================

for roi_idx = 1:numel(roi_list)

    roiName = roi_list{roi_idx};

    fprintf('\n');
    fprintf('############################################################\n');
    fprintf('ROI: %s\n', roiName);
    fprintf('############################################################\n');


    %% --------------------------------------------------------------------
    % Check ROI
    %% --------------------------------------------------------------------

    if ~isfield(AllTestMaps, roiName)
        warning('ROI %s is not found in AllTestMaps.', roiName);
        continue;
    end

    roiResults = AllTestMaps.(roiName);
    nSub = numel(roiResults);


    %% --------------------------------------------------------------------
    % Subject x method x component array
    %% --------------------------------------------------------------------

    absCorr = nan( ...
        nSub, ...
        numel(fieldNames), ...
        Kcomp);


    %% ====================================================================
    % Compute subject-level component-wise absolute correlations
    %% ====================================================================

    for s = 1:nSub

        subj = roiResults(s).subject;

        fprintf('Processing %s | ROI = %s\n', subj, roiName);

        seed = roiResults(s).rseedtest_vals;


        for m = 1:numel(fieldNames)

            fieldName = fieldNames{m};

            if ~isfield(roiResults(s), fieldName)
                error( ...
                    'Field "%s" is missing for %s | %s.', ...
                    fieldName, subj, roiName);
            end

            maps = roiResults(s).(fieldName);

            if size(maps, 2) < Kcomp
                error( ...
                    ['Only %d components are available for %s | %s | %s, ' ...
                     'but Kcomp = %d.'], ...
                    size(maps,2), ...
                    subj, ...
                    roiName, ...
                    methodNames{m}, ...
                    Kcomp);
            end


            for k = 1:Kcomp

                r = corr( ...
                    seed, ...
                    maps(:,k), ...
                    'Type', 'Pearson');

                absCorr(s,m,k) = abs(r);

            end

        end


        %% ----------------------------------------------------------------
        % Store subject-level correlations
        %% ----------------------------------------------------------------

        currentCorr = squeeze(absCorr(s,:,:));   % method x component

        ComponentWiseCorrResults.(roiName).subject(s).name = ...
            subj;

        ComponentWiseCorrResults.(roiName).subject(s).absCorr = ...
            currentCorr;

        ComponentWiseCorrResults.(roiName).subject(s).table = ...
            array2table( ...
                currentCorr, ...
                'VariableNames', {'K1','K2','K3','K4'}, ...
                'RowNames', methodNames);

    end


    %% ====================================================================
    % Group mean and standard error across subjects
    %% ====================================================================

    % Dimensions after squeeze:
    %   method x component

    meanAbsCorr = squeeze( ...
        mean(absCorr, 1, 'omitnan'));

    stdAbsCorr = squeeze( ...
        std(absCorr, 0, 1, 'omitnan'));

    nValid = squeeze( ...
        sum(~isnan(absCorr), 1));

    seAbsCorr = ...
        stdAbsCorr ./ sqrt(nValid);


    %% ====================================================================
    % Store ROI-level numerical results
    %% ====================================================================

    ComponentWiseCorrResults.(roiName).absCorr = ...
        absCorr;

    ComponentWiseCorrResults.(roiName).subjects = ...
        {roiResults.subject};

    ComponentWiseCorrResults.(roiName).methodNames = ...
        methodNames;

    ComponentWiseCorrResults.(roiName).Kcomp = ...
        Kcomp;

    ComponentWiseCorrResults.(roiName).meanAbsCorr = ...
        meanAbsCorr;

    ComponentWiseCorrResults.(roiName).stdAbsCorr = ...
        stdAbsCorr;

    ComponentWiseCorrResults.(roiName).seAbsCorr = ...
        seAbsCorr;


    %% ====================================================================
    % Display group-level numerical results
    %% ====================================================================

    fprintf('\nMean absolute correlations | ROI = %s\n', roiName);

    Tmean = array2table( ...
        meanAbsCorr, ...
        'VariableNames', {'K1','K2','K3','K4'}, ...
        'RowNames', methodNames);

    disp(Tmean);

    fprintf('Standard errors | ROI = %s\n', roiName);

    Tse = array2table( ...
        seAbsCorr, ...
        'VariableNames', {'K1','K2','K3','K4'}, ...
        'RowNames', methodNames);

    disp(Tse);


    %% ====================================================================
    % Plot group mean +/- SE
    %% ====================================================================

    fig = figure( ...
        'Name', sprintf( ...
            '%s_ComponentWise_HeldOutCorrelation', ...
            roiName), ...
        'Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1 1 8.5 4.5]);

    hold on;


    for m = 1:numel(methodNames)

        xCurrent = xBase + xShift(m);

        errorbar( ...
            xCurrent, ...
            meanAbsCorr(m,:), ...
            seAbsCorr(m,:), ...
            'LineStyle', 'none', ...
            'Marker', markerList{m}, ...
            'MarkerSize', 12, ...
            'LineWidth', 1.2, ...
            'CapSize', 8);

    end


    %% --------------------------------------------------------------------
    % Figure formatting
    %% --------------------------------------------------------------------

    xlabel( ...
        'Canonical component', ...
        'FontSize', 22);

    ylabel( ...
        'Absolute correlation', ...
        'FontSize', 22);

    xticks(1:Kcomp);
    xticklabels({'1','2','3','4'});

    xlim([0.7 4.3]);
    ylim([0 0.65]);

    legend( ...
        methodNames, ...
        'Location', 'northoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off');

    box on;
    grid off;

    set( ...
        gca, ...
        'FontSize', 20, ...
        'LineWidth', 1);


    %% ====================================================================
    % Save figure
    %% ====================================================================

    pdfFile = fullfile( ...
        figureDir, ...
        sprintf( ...
            '%s_ComponentWise_HeldOutSeedCorrelation.pdf', ...
            roiName));

    pngFile = fullfile( ...
        figureDir, ...
        sprintf( ...
            '%s_ComponentWise_HeldOutSeedCorrelation.png', ...
            roiName));

    exportgraphics( ...
        fig, ...
        pdfFile, ...
        'ContentType', 'vector');

    exportgraphics( ...
        fig, ...
        pngFile, ...
        'Resolution', 300);

    fprintf('Figures saved:\n');
    fprintf('  %s\n', pdfFile);
    fprintf('  %s\n', pngFile);

end


%% ========================================================================
% Save numerical results
% =========================================================================

outputFile = fullfile( ...
    outputDir, ...
    'ComponentWiseCorrResults.mat');

save( ...
    outputFile, ...
    'ComponentWiseCorrResults', ...
    'roi_list', ...
    'methodNames', ...
    'fieldNames', ...
    'Kcomp', ...
    '-v7.3');


%% ========================================================================
% Completion message
% =========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Component-wise held-out correlation analysis completed.\n');
fprintf('============================================================\n');

fprintf('\nNumerical results saved to:\n%s\n', outputFile);
fprintf('\nFigures saved to:\n%s\n', figureDir);
fprintf('\nDone.\n');
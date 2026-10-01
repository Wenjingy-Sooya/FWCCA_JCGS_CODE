%% ============================================================
% z0 Sensitivity Analysis:
% Oracle-referenced canonical-angle evaluation
%
% Purpose
% -------
% This script evaluates the sensitivity of the FWCCA canonical directions
% to different initial activation indicators z0 by computing canonical
% angles relative to the Oracle initialization.
%
% Input
% -----
%   Results/FWCCA_z0Sensitivity_AllSNR.mat
%
% Required structure
% ------------------
%   z0SensitivityResults(sidx).SNR_dB
%   z0SensitivityResults(sidx).fixed_h
%   z0SensitivityResults(sidx).GFCCAresults_currentSNR
%   z0SensitivityResults(sidx).LGFCCAresults_currentSNR
%   z0SensitivityResults(sidx).GLFCCAresults_currentSNR
%
% The z0 initialization ordering in each ResultCell is assumed to be:
%
%   1. Oracle
%   2. Noisy oracle
%   3. Preliminary CCA
%   4. Random matched
%
% Evaluation
% ----------
% For each:
%
%   - SNR level
%   - fixed bandwidth h
%   - FWCCA method
%   - canonical component
%   - z0 initialization comparison
%
% the canonical angle relative to the Oracle initialization is computed
% across the Monte Carlo replicates and summarized using:
%
%   - mean
%   - median
%   - standard deviation
%   - 95th percentile
%   - 99th percentile
%   - maximum
%
% Output
% ------
%   Results/
%       FWCCA_z0Sensitivity_oracleReferenced_angle_summary_all.mat
%
%   Results/
%       FWCCA_z0Sensitivity_oracleReferenced_angle_summary_all.csv
%
%% ============================================================
clear; clc;
%% ============================================================
% Project directories
%% ============================================================

% Directory containing this script:
%   fMRISimulation/02_z0SensitivityAnalysis/
script_dir = fileparts(mfilename('fullpath'));

% Project root directory:
%   fMRISimulation/
project_root = fileparts(script_dir);

% Shared FWCCA functions:
%   fMRISimulation/Functions/
functions_dir = fullfile(project_root, 'Functions');

if ~exist(functions_dir, 'dir')
    error('Functions directory does not exist:\n%s', functions_dir);
end

% Add the project root so that the MATLAB package
%   fMRISimulation/+z0sensitivity/
% can be accessed as z0sensitivity.*
addpath(project_root);

% Add shared FWCCA functions
addpath(functions_dir);


%% ============================================================
% Results directory
%% ============================================================

% Input and output results:
%   fMRISimulation/02_z0SensitivityAnalysis/Results/
results_dir = fullfile(script_dir, 'Results');

if ~exist(results_dir, 'dir')
    error('Results directory does not exist:\n%s', results_dir);
end


%% ============================================================
% Load z0 sensitivity simulation results
%% ============================================================

input_file = fullfile( ...
    results_dir, ...
    'FWCCA_z0Sensitivity_AllSNR.mat');

if ~exist(input_file, 'file')
    error('Simulation result file does not exist:\n%s', input_file);
end

load(input_file, 'z0SensitivityResults');

fprintf('\n============================================================\n');
fprintf('Loaded z0 sensitivity simulation results from:\n%s\n', input_file);
fprintf('Number of SNR levels: %d\n', numel(z0SensitivityResults));
fprintf('============================================================\n');

%% ============================================================
% output directories
%% ============================================================
output_dir = fullfile(script_dir, 'Results');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


%% ============================================================
% Settings
%% ============================================================

methodNames = { ...
    'G-FWCCA', ...
    'LG-FWCCA', ...
    'GL-FWCCA'};

resultVarNames = { ...
    'GFCCAresults_currentSNR', ...
    'LGFCCAresults_currentSNR', ...
    'GLFCCAresults_currentSNR'};

comparisonNames = { ...
    'Oracle vs noisy oracle', ...
    'Oracle vs preliminary CCA', ...
    'Oracle vs random matched'};

nMethods = numel(methodNames);
nComparisons = numel(comparisonNames);

% Number of canonical components used in the simulation
Kcomp = 2;

% Container for summary rows
rows = {};


%% ============================================================
% Loop over SNR levels
%% ============================================================

for sidx = 1:numel(z0SensitivityResults)

    %% --------------------------------------------------------
    % Current SNR-specific result structure
    %% --------------------------------------------------------

    S = z0SensitivityResults(sidx);

    SNR_dB = S.SNR_dB;
    fixed_h = S.fixed_h;

    fprintf('\n----------------------------------------\n');
    fprintf('SNR = %+d dB | fixed h = %.2f\n', ...
        SNR_dB, fixed_h);
    fprintf('----------------------------------------\n');


    %% ========================================================
    % Loop over FWCCA methods
    %% ========================================================

    for methodIdx = 1:nMethods

        resultVarName = resultVarNames{methodIdx};


        %% ----------------------------------------------------
        % Check that result field exists
        %% ----------------------------------------------------

        if ~isfield(S, resultVarName)

            error( ...
                'Cannot find field "%s" for SNR = %+d dB.', ...
                resultVarName, SNR_dB);

        end


        %% ----------------------------------------------------
        % Extract current FWCCA results
        %
        % Expected dimensions:
        %
        %   M x 4
        %
        % Columns correspond to:
        %   1 Oracle
        %   2 Noisy oracle
        %   3 Preliminary CCA
        %   4 Random matched
        %% ----------------------------------------------------

        ResultCell = S.(resultVarName);

        fprintf('  %s\n', methodNames{methodIdx});


        %% ----------------------------------------------------
        % Basic structure check
        %% ----------------------------------------------------

        if size(ResultCell,2) ~= 4

            error( ...
                ['Expected 4 z0 initialization columns for ' ...
                 '%s at SNR = %+d dB, but found %d.'], ...
                methodNames{methodIdx}, ...
                SNR_dB, ...
                size(ResultCell,2));

        end


        %% ----------------------------------------------------
        % Compute oracle-referenced canonical angles
        %
        % angles(m,c,k)
        %
        % m = Monte Carlo run
        %
        % c = comparison:
        %       1 Oracle vs noisy oracle
        %       2 Oracle vs preliminary CCA
        %       3 Oracle vs random matched
        %
        % k = canonical component
        %% ----------------------------------------------------

        angles = z0sensitivity.computeCanonicalAngles( ...
            ResultCell, ...
            Kcomp);


        %% ====================================================
        % Summarize each component and comparison separately
        %% ====================================================

        for k = 1:Kcomp

            for c = 1:nComparisons

                %% --------------------------------------------
                % One angle per Monte Carlo replicate
                %% --------------------------------------------

                angleVals = angles(:,c,k);


                %% --------------------------------------------
                % Remove NaN
                %
                % Non-converged or unavailable pairs are
                % represented by NaN.
                %% --------------------------------------------

                angleVals = angleVals(~isnan(angleVals));

                nAngles = numel(angleVals);


                %% --------------------------------------------
                % Summary statistics
                %% --------------------------------------------

                if nAngles == 0

                    meanAngle   = NaN;
                    medianAngle = NaN;
                    sdAngle     = NaN;

                    p95Angle = NaN;
                    p99Angle = NaN;

                    maxAngle = NaN;

                else

                    meanAngle = mean(angleVals);

                    medianAngle = median(angleVals);

                    sdAngle = std(angleVals);

                    p95Angle = prctile( ...
                        angleVals, ...
                        95);

                    p99Angle = prctile( ...
                        angleVals, ...
                        99);

                    maxAngle = max(angleVals);

                end


                %% --------------------------------------------
                % Store one summary row
                %% --------------------------------------------

                rows(end+1,:) = { ...
                    SNR_dB, ...
                    fixed_h, ...
                    methodNames{methodIdx}, ...
                    k, ...
                    comparisonNames{c}, ...
                    meanAngle, ...
                    medianAngle, ...
                    sdAngle, ...
                    p95Angle, ...
                    p99Angle, ...
                    maxAngle, ...
                    nAngles};

            end

        end

    end

end


%% ============================================================
% Create summary table
%% ============================================================

AngleComparisonTable = cell2table( ...
    rows, ...
    'VariableNames', { ...
        'SNR', ...
        'fixed_h', ...
        'Method', ...
        'Component', ...
        'Comparison', ...
        'meanAngle', ...
        'medianAngle', ...
        'sdAngle', ...
        'p95Angle', ...
        'p99Angle', ...
        'maxAngle', ...
        'nAngles'});


%% ============================================================
% Sort table
%% ============================================================

AngleComparisonTable = sortrows( ...
    AngleComparisonTable, ...
    {'SNR','Method','Component','Comparison'});


%% ============================================================
% Display summary table
%% ============================================================

disp(AngleComparisonTable);


%% ============================================================
% Basic checks
%% ============================================================

expectedRows = ...
    numel(z0SensitivityResults) * ...
    nMethods * ...
    Kcomp * ...
    nComparisons;

fprintf('\nNumber of summary rows: %d\n', ...
    height(AngleComparisonTable));

fprintf('Expected number of rows: %d\n', ...
    expectedRows);


if height(AngleComparisonTable) ~= expectedRows

    warning( ...
        'Unexpected number of rows in AngleComparisonTable.');

end


%% ============================================================
% Check expected SNR levels and fixed bandwidths
%% ============================================================

SNR_check = [z0SensitivityResults.SNR_dB];
h_check = [z0SensitivityResults.fixed_h];

fprintf('\nSNR-dependent bandwidths:\n');

BandwidthTable = table( ...
    SNR_check(:), ...
    h_check(:), ...
    'VariableNames', ...
    {'SNR_dB','fixed_h'});

disp(BandwidthTable);


%% ============================================================
% Save outputs
%% ============================================================

outputMat = fullfile( ...
    output_dir, ...
    'FWCCA_z0Sensitivity_oracleReferenced_angle_summary_all.mat');

outputCSV = fullfile( ...
    output_dir, ...
    'FWCCA_z0Sensitivity_oracleReferenced_angle_summary_all.csv');


save( ...
    outputMat, ...
    'AngleComparisonTable');


writetable( ...
    AngleComparisonTable, ...
    outputCSV);


fprintf('\nSaved:\n');
fprintf('  %s\n', outputMat);
fprintf('  %s\n', outputCSV);

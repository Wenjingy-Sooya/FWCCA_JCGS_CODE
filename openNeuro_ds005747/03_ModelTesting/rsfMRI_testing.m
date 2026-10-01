%% ========================================================================
% rsfMRI_Testing.m
%
% Purpose:
%   Perform held-out testing for the resting-state fMRI FWCCA analysis.
%
%   For each subject and each default mode network (DMN) ROI:
%
%       AnG  : Angular Gyrus
%       PCgG : Posterior Cingulate Gyrus
%       PCu  : Precuneus
%       SFG  : Superior Frontal Gyrus
%
%   the script:
%
%   1. Loads the training and held-out testing inputs prepared in Step 01
%      (rsfMRI_FWCCA_Inputs.mat).
%
%   2. Loads the LG-FWCCA and GL-FWCCA hyperparameters selected by
%      cross-validation in Step 02.
%
%   3. Fits CCA and G-FWCCA to the complete training data, and refits
%      LG-FWCCA and GL-FWCCA using the CV-selected (h,r) values for
%      K = 1, 2, 3, and 4 canonical components.
%
%   4. Projects the independent held-out test data onto the canonical
%      directions estimated from the training data.
%
%   5. Uses the held-out seed-correlation map (rseedtest_vals), prepared
%      independently from the held-out data in Step 01, as the reference
%      map for testing.
%
%   6. Computes the absolute Pearson correlation between the held-out
%      seed-correlation map and each held-out canonical spatial map.
%
%   7. Computes the cumulative absolute correlation across the first
%      K = 1, 2, 3, and 4 canonical components.
%
% Inputs:
%   ../01_PrepareAnalysisInputs/Results/rsfMRI_FWCCA_Inputs.mat
%   ../02_HyperparameterSelection/Results/rsfMRI_LG_GL_CV_Results.mat
%
% Outputs:
%   Results/rsfMRI_TestingResults.mat
%
%   The saved results contain:
%
%   AllTestCorrResults
%       Subject-level cumulative held-out correlations for CCA,
%       G-FWCCA, LG-FWCCA, and GL-FWCCA.
%
%   AllTestMaps
%       Held-out seed-correlation maps and canonical spatial maps for
%       the four methods.
%
%   GroupTestMeans
%       ROI-level mean cumulative correlations across subjects.
%
%   A group-level summary table is also generated for subsequent
%   evaluation and visualization in Step 04.
%
% IMPORTANT:
%   The held-out test data are not used for model fitting or
%   hyperparameter selection. All canonical directions are estimated
%   using the training data only.
%
% =========================================================================
% =========================================================================

clear;
clc;
close all;


%% ========================================================================
% Project paths
% =========================================================================

% Directory containing this script:
%   openNeuro_ds005747/03_ModelTesting
scriptDir = fileparts(mfilename('fullpath'));

% Dataset project directory:
%   openNeuro_ds005747
datasetRoot = fileparts(scriptDir);

% Step 01 output directory:
%   openNeuro_ds005747/01_PrepareAnalysisInputs/Results
inputRoot = fullfile( ...
    datasetRoot, ...
    '01_PrepareAnalysisInputs', ...
    'Results');

% Step 02 output directory:
%   openNeuro_ds005747/02_HyperparameterSelection/Results
cvRoot = fullfile( ...
    datasetRoot, ...
    '02_HyperparameterSelection', ...
    'Results');

% Shared FWCCA functions:
%   openNeuro_ds005747/Functions
functionsDir = fullfile( ...
    datasetRoot, ...
    'Functions');


%% ========================================================================
% Check required directories
% =========================================================================

if ~exist(functionsDir, 'dir')
    error('Functions directory not found:\n%s', functionsDir);
end

addpath(functionsDir);


%% ========================================================================
% Load Step 01: prepared training and held-out testing inputs
% =========================================================================

inputFile = fullfile( ...
    inputRoot, ...
    'rsfMRI_FWCCA_Inputs.mat');

if ~isfile(inputFile)
    error(['Step 01 input file not found:\n%s\n\n' ...
           'Run Prepare_rsfMRI_FWCCA_Inputs.m first.'], ...
           inputFile);
end

fprintf('Loading Step 01 inputs:\n%s\n\n', inputFile);

S_input = load( ...
    inputFile, ...
    'rsfMRIInputs');

if ~isfield(S_input, 'rsfMRIInputs')
    error('Variable "rsfMRIInputs" was not found in:\n%s', inputFile);
end

rsfMRIInputs = S_input.rsfMRIInputs;


%% ========================================================================
% Load Step 02: LG-FWCCA / GL-FWCCA cross-validation results
% =========================================================================

cvFile = fullfile( ...
    cvRoot, ...
    'rsfMRI_LG_GL_CV_Results.mat');

if ~isfile(cvFile)
    error(['Step 02 cross-validation result file not found:\n%s\n\n' ...
           'Run LG_GL_FWCCA_CrossValidation.m first.'], ...
           cvFile);
end

fprintf('Loading Step 02 CV results:\n%s\n\n', cvFile);

S_cv = load( ...
    cvFile, ...
    'rsfMRI_CV_Results');

if ~isfield(S_cv, 'rsfMRI_CV_Results')
    error('Variable "rsfMRI_CV_Results" was not found in:\n%s', cvFile);
end

rsfMRICVResults = S_cv.rsfMRI_CV_Results;


%% ========================================================================
% Input loading completed
% =========================================================================

fprintf('============================================================\n');
fprintf('Step 01 inputs and Step 02 CV results loaded successfully.\n');
fprintf('============================================================\n\n');

%% ============================================================
% Analysis settings
%% ============================================================
Kmaxcomp = 4;
svdTol  = 1e-1;
ridgeTol= 1e-4;

%% ========================================================================
% Subjects and DMN ROIs
%% ============================================================
sub_list = {'sub-011', 'sub-012', 'sub-018',  'sub-021', 'sub-022'};
roi_list = {'AnG', 'PCgG', 'PCu',  'SFG'};

%% ============================================================
% Store testing results
%% ============================================================

AllTestCorrResults = struct();
AllTestMaps = struct();

%% ========================================================================
% Held-out testing: ROI x subject
% =========================================================================

for roi_idx = 1:numel(roi_list)

    roiName = roi_list{roi_idx};

    fprintf('\n############################################################\n');
    fprintf('ROI: %s\n', roiName);
    fprintf('############################################################\n');

    % Check that Step 01 inputs and Step 02 CV results are available
    if ~isfield(rsfMRIInputs, roiName)
        warning('ROI %s is not found in rsfMRIInputs.', roiName);
        continue;
    end

    if ~isfield(rsfMRICVResults, roiName)
        warning('ROI %s is not found in rsfMRICVResults.', roiName);
        continue;
    end

    Inputs    = rsfMRIInputs.(roiName);
    CVresults = rsfMRICVResults.(roiName);


    %% ====================================================================
    % Subject loop
    % =====================================================================

    for sub_id = 1:numel(sub_list)

        subj = sub_list{sub_id};

        fprintf('\n============================================================\n');
        fprintf('Testing %s | ROI = %s\n', subj, roiName);
        fprintf('============================================================\n');


        %% ----------------------------------------------------------------
        % Load stored training and held-out testing inputs
        % -----------------------------------------------------------------

        X1seed_train   = Inputs(sub_id).X1seed_train;
        X2seed_train   = Inputs(sub_id).X2seed_train;
        Xseed_test     = Inputs(sub_id).Xseed_test;
        Wg_0           = Inputs(sub_id).Wg_0;
        rseedtest_vals = Inputs(sub_id).rseedtest_vals;

        if ~isfield(Inputs(sub_id), 'TR')
            error(['TR is not stored in Inputs(%d) for %s | %s.\n' ...
                   'Please check fieldnames(Inputs(sub_id)).'], ...
                   sub_id, subj, roiName);
        end

        TR     = Inputs(sub_id).TR;
        Ttrain = size(X1seed_train, 2);
        I_T    = eye(Ttrain);


        %% ----------------------------------------------------------------
        % Fit CCA and G-FWCCA on the complete training data
        % -----------------------------------------------------------------

        resCCA = Run_FWCCA_Family( ...
            X1seed_train, X2seed_train, ...
            I_T, I_T, I_T, I_T, ...
            Kmaxcomp, 100, 1e-4, ...
            'prescale', ridgeTol, svdTol);

        resG = Run_FWCCA_Family( ...
            X1seed_train, X2seed_train, ...
            Wg_0, Wg_0, I_T, I_T, ...
            Kmaxcomp, 100, 1e-4, ...
            'prescale', ridgeTol, svdTol);


        %% ----------------------------------------------------------------
        % Project independent held-out test data
        % -----------------------------------------------------------------

        CCA_test = (Xseed_test - resCCA.X1_localMean) * resCCA.Ak_raw;
        G_test   = (Xseed_test - resG.X1_localMean)   * resG.Ak_raw;

        AllTestMaps.(roiName)(sub_id).subject        = subj;
        AllTestMaps.(roiName)(sub_id).rseedtest_vals = rseedtest_vals;
        AllTestMaps.(roiName)(sub_id).CCA            = CCA_test;
        AllTestMaps.(roiName)(sub_id).G_FWCCA        = G_test;


        %% ----------------------------------------------------------------
        % Component-wise held-out correlations for CCA and G-FWCCA
        % -----------------------------------------------------------------

        rseedtest_norm = normalize_to_unit_interval(rseedtest_vals);

        CCA_corr_Vec = arrayfun(@(j) corr( ...
            rseedtest_norm, ...
            normalize_to_unit_interval(CCA_test(:,j)), ...
            'Type', 'Pearson'), 1:Kmaxcomp);

        GCCA_corr_Vec = arrayfun(@(j) corr( ...
            rseedtest_norm, ...
            normalize_to_unit_interval(G_test(:,j)), ...
            'Type', 'Pearson'), 1:Kmaxcomp);


        %% ----------------------------------------------------------------
        % CV-selected LG-FWCCA and GL-FWCCA parameters
        %
        % Columns:
        %   1 = K, 2 = h, 3 = r
        % -----------------------------------------------------------------

        LG_BestParaTable = CVresults(sub_id).LG_BestParamTable;
        GL_BestParaTable = CVresults(sub_id).GL_BestParamTable;

        CCA_cumsum  = nan(1, Kmaxcomp);
        GCCA_cumsum = nan(1, Kmaxcomp);
        LG_cumsum   = nan(1, Kmaxcomp);
        GL_cumsum   = nan(1, Kmaxcomp);


        %% =================================================================
        % Refit LG-FWCCA and GL-FWCCA for K = 1,...,Kmaxcomp
        % ==================================================================

        for k = 1:Kmaxcomp

            fprintf('\n------------------------------------------------------------\n');
            fprintf('%s | %s | K = %d\n', subj, roiName, k);
            fprintf('------------------------------------------------------------\n');

            % CV-selected parameters for the current K
            LG_h = LG_BestParaTable{k, 2};
            LG_r = LG_BestParaTable{k, 3};
            GL_h = GL_BestParaTable{k, 2};
            GL_r = GL_BestParaTable{k, 3};

            fprintf('LG selected: h = %.4f | r = %.4f\n', LG_h, LG_r);
            fprintf('GL selected: h = %.4f | r = %.4f\n', GL_h, GL_r);


            %% -------------------------------------------------------------
            % Construct temporal local kernels
            % --------------------------------------------------------------

            K_LG = makeLocalKernel( ...
                Ttrain, LG_r, 'gaussian', LG_h, TR);

            K_GL = makeLocalKernel( ...
                Ttrain, GL_r, 'gaussian', GL_h, TR);


            %% -------------------------------------------------------------
            % Refit LG-FWCCA and GL-FWCCA on the complete training data
            % --------------------------------------------------------------

            resLG = Run_FWCCA_Family( ...
                X1seed_train, X2seed_train, ...
                Wg_0, Wg_0, K_LG, K_LG, ...
                k, 100, 1e-4, ...
                'postscale', ridgeTol, svdTol);

            resGL = Run_FWCCA_Family( ...
                X1seed_train, X2seed_train, ...
                Wg_0, Wg_0, K_GL, K_GL, ...
                k, 100, 1e-4, ...
                'prescale', ridgeTol, svdTol);


            %% -------------------------------------------------------------
            % Project independent held-out test data
            % --------------------------------------------------------------

            LG_test = (Xseed_test - resLG.X1_localMean) * resLG.Ak_raw;
            GL_test = (Xseed_test - resGL.X1_localMean) * resGL.Ak_raw;


            %% -------------------------------------------------------------
            % Component-wise held-out correlations
            % --------------------------------------------------------------

            LGcca_corr_Vec = arrayfun(@(j) corr( ...
                rseedtest_norm, ...
                normalize_to_unit_interval(LG_test(:,j)), ...
                'Type', 'Pearson'), 1:k);

            GLcca_corr_Vec = arrayfun(@(j) corr( ...
                rseedtest_norm, ...
                normalize_to_unit_interval(GL_test(:,j)), ...
                'Type', 'Pearson'), 1:k);


            %% -------------------------------------------------------------
            % Cumulative absolute correlations
            % --------------------------------------------------------------

            CCA_cumsum(k)  = sum(abs(CCA_corr_Vec(1:k)),  'omitnan');
            GCCA_cumsum(k) = sum(abs(GCCA_corr_Vec(1:k)), 'omitnan');
            LG_cumsum(k)   = sum(abs(LGcca_corr_Vec),      'omitnan');
            GL_cumsum(k)   = sum(abs(GLcca_corr_Vec),      'omitnan');


            % Store Kmax canonical maps for subsequent spatial evaluation
            if k == Kmaxcomp
                AllTestMaps.(roiName)(sub_id).LG_FWCCA = LG_test;
                AllTestMaps.(roiName)(sub_id).GL_FWCCA = GL_test;
            end

        end


        %% ----------------------------------------------------------------
        % Store subject-level held-out testing results
        % -----------------------------------------------------------------

        AllTestCorrResults.(roiName)(sub_id).subject     = subj;
        AllTestCorrResults.(roiName)(sub_id).CCA_cumsum  = CCA_cumsum;
        AllTestCorrResults.(roiName)(sub_id).GCCA_cumsum = GCCA_cumsum;
        AllTestCorrResults.(roiName)(sub_id).LG_cumsum   = LG_cumsum;
        AllTestCorrResults.(roiName)(sub_id).GL_cumsum   = GL_cumsum;

    end

end

%% ========================================================================
% Compute group means across subjects
% =========================================================================

GroupTestMeans = struct();

methods = {'CCA_cumsum', 'GCCA_cumsum', 'LG_cumsum', 'GL_cumsum'};

for roi_idx = 1:numel(roi_list)

    roiName = roi_list{roi_idx};

    if ~isfield(AllTestCorrResults, roiName)
        continue;
    end

    groupResults = AllTestCorrResults.(roiName);

    for m = 1:numel(methods)

        methodName = methods{m};

        % Subjects x K cumulative-correlation matrix
        mat = cell2mat(arrayfun( ...
            @(s) s.(methodName)(:)', ...
            groupResults(:), ...
            'UniformOutput', false));

        % Group mean for each cumulative K
        GroupTestMeans.(roiName).(methodName) = ...
            mean(mat, 1, 'omitnan');

    end

end


%% ========================================================================
% Display group means
% =========================================================================

for roi_idx = 1:numel(roi_list)

    roiName = roi_list{roi_idx};

    if ~isfield(GroupTestMeans, roiName)
        continue;
    end

    fprintf('\n############################################################\n');
    fprintf('Group mean | ROI = %s\n', roiName);
    fprintf('############################################################\n');

    fprintf('CCA:\n');
    disp(GroupTestMeans.(roiName).CCA_cumsum);

    fprintf('G-FWCCA:\n');
    disp(GroupTestMeans.(roiName).GCCA_cumsum);

    fprintf('LG-FWCCA:\n');
    disp(GroupTestMeans.(roiName).LG_cumsum);

    fprintf('GL-FWCCA:\n');
    disp(GroupTestMeans.(roiName).GL_cumsum);

end


%% ============================================================
% Save testing results
%% ============================================================

output_dir = fullfile(scriptDir, 'Results');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end



outputFile = fullfile( ...
    output_dir, ...
   'rsfMRI_TestingResults.mat');

save( ...
    outputFile, ...
    'AllTestCorrResults', ...
    'GroupTestMeans', ...
    'AllTestMaps', ...
    'roi_list', ...
    'sub_list', ...
    'Kmaxcomp', ...
    'svdTol', ...
    '-v7.3');

fprintf('\nTesting results saved to:\n%s\n', outputFile);

%% ============================================================
% Save group-level summary table
%% ============================================================


summaryRows = {};

for roi_idx = 1:numel(roi_list)

    roiName = roi_list{roi_idx};

    if ~isfield(GroupTestMeans, roiName)
        continue;
    end
 
    CCA_mean =  GroupTestMeans.(roiName).CCA_cumsum;
    GCCA_mean =  GroupTestMeans.(roiName).GCCA_cumsum;
    LG_mean = GroupTestMeans.(roiName).LG_cumsum;
    GL_mean = GroupTestMeans.(roiName).GL_cumsum;

    for k = 1:Kmaxcomp

        summaryRows(end+1, :) = { ...
            roiName, ...
            k, ...
            CCA_mean(k),...
            GCCA_mean(k),...
            LG_mean(k), ...
            GL_mean(k)};

    end

end

GroupSummaryTable = cell2table( ...
    summaryRows, ...
    'VariableNames', { ...
        'ROI', ...
        'K', ...
        'CCA',...
        'G-FWCCA',...
        'LG_FWCCA', ...
        'GL_FWCCA'});

excelFile = fullfile( ...
    output_dir, ...
    'rsfMRI_TestingGroupSummary_%s.xlsx');

writetable(GroupSummaryTable, excelFile);

fprintf('Group summary Excel saved to:\n%s\n', excelFile);


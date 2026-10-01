% =========================================================================
% LG_GL_FWCCA_CrossValidation.m
%
% Purpose:
%   Select the local-kernel hyperparameters (h,r) for LG-FWCCA and
%   GL-FWCCA using the training data prepared in Step 01.
%
%   Cross-validation is performed separately for each subject and each
%   default mode network (DMN) ROI:
%
%       AnG  : Angular Gyrus
%       PCgG : Posterior Cingulate Gyrus
%       PCu  : Precuneus
%       SFG  : Superior Frontal Gyrus
%
%   Candidate temporal bandwidths h and spatial neighborhood radii r are
%   evaluated using voxel-wise cross-validation. Hyperparameter selection
%   is based on held-out seed-map correlation within the training data.
%
%   Only the training data are used in this script. The held-out testing
%   data stored in rsfMRI_FWCCA_Inputs.mat are not accessed.
%
% Input:
%   ../01_PrepareAnalysisInputs/Results/rsfMRI_FWCCA_Inputs.mat
%
% Output:
%   Results/rsfMRI_LG_GL_CV_Results.mat
%
% =========================================================================

clear;
clc;
close all;


%% ========================================================================
% Project paths
% =========================================================================

% Directory containing this script:
%   openNeuro_ds005747/02_ModelTraining
scriptDir = fileparts(mfilename('fullpath'));

% Dataset project directory:
%   openNeuro_ds005747
datasetRoot = fileparts(scriptDir);

% Parent directory containing both openNeuro_ds005747 and spm_25
projectRoot = fileparts(datasetRoot);

% Shared FWCCA functions:
%   openNeuro_ds005747/Functions
functionsDir = fullfile(datasetRoot, 'Functions');

% SPM25 directory
spmDir = fullfile(projectRoot, 'spm_25');

% Input generated in Step 01
inputFile = fullfile( ...
    datasetRoot, ...
    '01_PrepareAnalysisInputs', ...
    'Results', ...
    'rsfMRI_FWCCA_Inputs.mat');

% Output directory
outputDir = fullfile(scriptDir, 'Results');


%% ========================================================================
% Check required paths
% =========================================================================

if ~exist(functionsDir, 'dir')
    error('Functions directory not found:\n%s', functionsDir);
end

if ~exist(spmDir, 'dir')
    error('SPM25 directory not found:\n%s', spmDir);
end

if ~isfile(inputFile)
    error(['Step 01 input file not found:\n%s\n\n' ...
           'Run Prepare_rsfMRI_FWCCA_Inputs.m first.'], ...
           inputFile);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end


%% ========================================================================
% Add required paths
% =========================================================================

addpath(functionsDir);
addpath(spmDir)
%% ========================================================================
% Load prepared analysis inputs
% =========================================================================

fprintf('============================================================\n');
fprintf('LG-FWCCA / GL-FWCCA Cross-Validation\n');
fprintf('============================================================\n\n');

fprintf('Loading Step 01 inputs:\n%s\n\n', inputFile);

S = load(inputFile, 'rsfMRIInputs');

if ~isfield(S, 'rsfMRIInputs')
    error('Variable "rsfMRIInputs" was not found in:\n%s', inputFile);
end

rsfMRIInputs = S.rsfMRIInputs;


%% ========================================================================
% Neuromorphometrics atlas
% =========================================================================

atlas_file = fullfile( ...
    spmDir, ...
    'atlas', ...
    'Neuromorphometrics', ...
    'Neuromorphometrics.nii');

if ~isfile(atlas_file)
    error('Neuromorphometrics atlas not found:\n%s', atlas_file);
end

atlas_vol  = spm_vol(atlas_file);
atlas_data = spm_read_vols(atlas_vol);


%% ========================================================================
% Cross-validation settings
% =========================================================================
% Maximum number of canonical components
Kmaxcomp = 4;
% Candidate spatial neighborhood radii
r_list = 0:1:5;
% Candidate temporal-kernel bandwidths
h_list = 1:0.5:5;
% Number of voxel-wise cross-validation folds
num_fold = 10;
% Ridge and SVD tolerance used in FWCCA fitting
svdTol  = 1e-1;
ridgeTol= 1e-4;

%% ========================================================================
% Subjects and DMN ROIs
%% ============================================================
sub_list = {'sub-011', 'sub-012', 'sub-018',  'sub-021', 'sub-022'};
group_name = {'AnG', 'PCgG', 'PCu',  'SFG'};
%% ============================================================



%% ========================================================================
% Initialise output structure
% =========================================================================

rsfMRI_CV_Results = struct();


%% ============================================================
% Main loop: ROI x subject
%% ============================================================

for id_l = 1:numel(group_name)

    gname = group_name{id_l};

    fprintf('\n');
    fprintf('############################################################\n');
    fprintf('ROI: %s\n', gname);
    fprintf('############################################################\n');


    %% --------------------------------------------------------
    % Retrieve stored inputs for current ROI
    %% --------------------------------------------------------

    if iscell(rsfMRIInputs)

        Inputs = rsfMRIInputs{id_l};

    elseif isstruct(rsfMRIInputs) && isfield(rsfMRIInputs,gname)

        Inputs = rsfMRIInputs.(gname);

    else

        error( ...
            'Cannot determine rsfMRIInputs structure for ROI %s.', ...
            gname);

    end

    ResultsPerGroup = struct([]);


    %% ========================================================
    % Subject loop
    %% ========================================================

    for sub_id = 1:numel(sub_list)

        subj = sub_list{sub_id};
        fprintf('\n');
        fprintf('============================================================\n');
        fprintf('Processing %s | ROI = %s\n', subj, gname);
        fprintf('============================================================\n');
       %% ====================================================
        % Recover stored inputs
        %% ====================================================
        X1seed_train = Inputs(sub_id).X1seed_train;
        X2seed_train = Inputs(sub_id).X2seed_train;
        Ttrain=size(X1seed_train,2);
        Wg_0 = Inputs(sub_id).Wg_0;
        seed_labels = Inputs(sub_id).seed_labels;
        TR = Inputs(sub_id).TR;    
        %% ====================================================
        % Reconstruct ROI voxel indices from atlas
        %% ====================================================
        seed_mask = ismember(atlas_data, seed_labels);
        seed_indices = find(seed_mask);
        %% ----------------------------------------------------
        % Check spatial consistency
        %% ----------------------------------------------------

        if numel(seed_indices) ~= size(X1seed_train,1)

            error(['ROI voxel-count mismatch for %s | %s:\n' ...
                 'Atlas ROI voxels       = %d\n' ...
                 'Rows in X1seed_train   = %d'], ...
                subj, ...
                gname, ...
                numel(seed_indices), ...
                size(X1seed_train,1));

        end


        fprintf('ROI voxels       = %d\n', numel(seed_indices));
        fprintf('Training points  = %d\n', Ttrain);
        fprintf('Kmaxcomp         = %d\n', Kmaxcomp);
        fprintf('TR               = %.3f sec\n', TR);
        fprintf('==== run local FWCCA varians  ====\n');

        %% ====================================================
        % LG-FWCCA: cross-validation for h and r
        %% ====================================================
      
        [LG_meanCumCorr, LG_BestParamTable, LG_CV_seedCorr] = ...
            CV_FWCCA_SeedCorrelation( ...
                seed_indices, ...
                X1seed_train, ...
                X2seed_train, ...
                Wg_0, ...
                Kmaxcomp, ...
                h_list, ...
                r_list, ...
                num_fold, ...
                TR, ...
                'postscale', ridgeTol, svdTol);

        fprintf('==== LG-FWCCA CV completed ====\n');
        disp(LG_BestParamTable);

        %% ====================================================
        % GL-FWCCA:
        % cross-validation for h and r
        %% ====================================================
        
        [GL_meanCumCorr, ...
         GL_BestParamTable, ...
         GL_CV_seedCorr] = ...
            CV_FWCCA_SeedCorrelation( ...
                seed_indices, ...
                X1seed_train, ...
                X2seed_train, ...
                Wg_0, ...
                Kmaxcomp, ...
                h_list, ...
                r_list, ...
                num_fold, ...
                TR, ...
                'prescale', ridgeTol, svdTol);

        fprintf('==== GL-FWCCA CV completed ====\n');
        disp(GL_BestParamTable);
        %% ====================================================
        % Store results
        %% ====================================================
        ResultsPerGroup(sub_id).subject =subj;
        ResultsPerGroup(sub_id).ROI = gname;
        %% ----------------------------------------------------
        % Stored analysis inputs
        %% ----------------------------------------------------

        ResultsPerGroup(sub_id).seed_labels = seed_labels;

        ResultsPerGroup(sub_id).seed_indices = ...
            seed_indices;

        ResultsPerGroup(sub_id).Wg_0 = ...
            Wg_0;

        ResultsPerGroup(sub_id).TR = ...
            TR;
        %% ----------------------------------------------------
        % Analysis settings
        %% ----------------------------------------------------

        ResultsPerGroup(sub_id).Kmaxcomp = ...
            Kmaxcomp;

        ResultsPerGroup(sub_id).h_list = ...
            h_list;

        ResultsPerGroup(sub_id).r_list = ...
            r_list;

        ResultsPerGroup(sub_id).num_fold = ...
            num_fold;


        %% ----------------------------------------------------
        % LG-FWCCA CV results
        %% ----------------------------------------------------

        ResultsPerGroup(sub_id).LG_meanCumCorr = ...
            LG_meanCumCorr;

        ResultsPerGroup(sub_id).LG_BestParamTable = ...
            LG_BestParamTable;

        ResultsPerGroup(sub_id).LG_CV_seedCorr = ...
            LG_CV_seedCorr;


        %% ----------------------------------------------------
        % GL-FWCCA CV results
        %% ----------------------------------------------------

        ResultsPerGroup(sub_id).GL_meanCumCorr = ...
            GL_meanCumCorr;

        ResultsPerGroup(sub_id).GL_BestParamTable = ...
            GL_BestParamTable;

        ResultsPerGroup(sub_id).GL_CV_seedCorr = ...
            GL_CV_seedCorr;

        ResultsPerGroup(sub_id).svdTol = svdTol;

        fprintf('\n');
        fprintf('%s | %s completed and stored.\n', subj,  gname);

    end

    %% --------------------------------------------------------
    % Store ROI-level results
    %% --------------------------------------------------------

    rsfMRI_CV_Results.(gname) = ResultsPerGroup;

end


%% ========================================================================
% Save cross-validation results
% =========================================================================

outputFile = fullfile( ...
    outputDir, ...
    'rsfMRI_LG_GL_CV_Results.mat');

save(outputFile, ...
    'rsfMRI_CV_Results', ...
    'r_list', ...
    'h_list', ...
    'Kmaxcomp', ...
    'num_fold', ...
    'svdTol', ...
    '-v7.3');


%% ========================================================================
% Completion message
% =========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('LG-FWCCA / GL-FWCCA cross-validation completed.\n');
fprintf('============================================================\n');

fprintf('\nResults saved to:\n%s\n', outputFile);
fprintf('\nDone.\n');


% =========================================================================
% Prepare_rsfMRI_FWCCA_Inputs.m
%
% Purpose:
%   Prepare the subject-level inputs used in the FWCCA analysis of the
%   resting-state fMRI dataset OpenNeuro ds005747.
%
%   The resting-state fMRI images used here are the preprocessed images
%   distributed with the dataset:
%
%       <subject>_task-rest_space-MNI305_preproc.nii
%
%   For each selected subject, the preprocessed fMRI image is first
%   resliced to the voxel grid of the Neuromorphometrics atlas. Four
%   atlas-defined default mode network (DMN) regions are then considered:
%
%       AnG  : Angular Gyrus
%       PCgG : Posterior Cingulate Gyrus
%       PCu  : Precuneus
%       SFG  : Superior Frontal Gyrus
%
%   For each ROI and subject, this script constructs the inputs required
%   by the subsequent CCA and FWCCA analyses.
%
%   No CCA/FWCCA model fitting or hyperparameter selection is performed
%   in this script.
%
% Output:
%   Results/rsfMRIInputs.mat
%
%   The output variable "rsfMRIInputs" contains the following fields for
%   each ROI and subject:
%
%       subject
%       X1seed_train
%       X2seed_train
%       Wg_0
%       Xseed_test
%       Kcomp
%       TR
%       seed_labels
%       rseedtrain_vals
%       rseedtest_vals
%
% =========================================================================

clear;
close all;
clc;

%% ========================================================================
% Project paths
% =========================================================================

% Directory containing this script:
%   openNeuro_ds005747/01_PrepareAnalysisInputs
scriptDir = fileparts(mfilename('fullpath'));

% Dataset project directory:
%   openNeuro_ds005747
datasetRoot = fileparts(scriptDir);

% Parent directory containing both openNeuro_ds005747 and spm_25
projectRoot = fileparts(datasetRoot);

% Dataset directory containing the subject folders:
%   openNeuro_ds005747/Data
dataRoot = fullfile(datasetRoot, 'Data');

% Shared project functions:
%   openNeuro_ds005747/Functions
functionsDir = fullfile(datasetRoot, 'Functions');

% SPM25 directory:
%   FWCCA_JCGS_CODE/spm_25
spmDir = fullfile(projectRoot, 'spm_25');


%% ========================================================================
% Check required directories
% =========================================================================

if ~exist(dataRoot, 'dir')
    error('Data directory not found:\n%s', dataRoot);
end

if ~exist(functionsDir, 'dir')
    error('Functions directory not found:\n%s', functionsDir);
end

if ~exist(spmDir, 'dir')
    error('SPM25 directory not found:\n%s', spmDir);
end


%% ========================================================================
% Add required paths
% =========================================================================

addpath(functionsDir);
addpath(spmDir);

spm('defaults', 'fmri');
spm_jobman('initcfg');

%% ========================================================================
% Output directory
% =========================================================================

outputDir = fullfile(scriptDir, 'Results');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

outputFile = fullfile(outputDir, 'rsfMRIInputs.mat');

%% ========================================================================
% Subjects included in the analysis
% =========================================================================

sub_list = {  ...
    'sub-011',...
    'sub-012',...
    'sub-018',...
    'sub-021', ...
    'sub-022'};



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
vol_size_3d = atlas_vol.dim;


%% ========================================================================
% DMN ROI definitions
% =========================================================================
%
% Left/right atlas labels are combined for each ROI.
%

group_labels = { ...
    [106 107], ...     % AnG
    [166 167], ...     % PCgG
    [168 169], ...     % PCu
    [190 191]};        % SFG

group_name = { ...
    'AnG', ...
    'PCgG', ...
    'PCu', ...
    'SFG'};


%% ========================================================================
% Analysis settings
% =========================================================================

% Maximum number of canonical components used in subsequent analyses
Kcomp = 4;


%% ========================================================================
% Prepare subject-level rs-fMRI inputs
% =========================================================================

rsfMRIInputs = struct();


for id_l = 1:numel(group_name)

    gname = group_name{id_l};
    seed_labels = group_labels{id_l};

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('ROI: %s\n', gname);
    fprintf('Atlas labels: %s\n', mat2str(seed_labels));
    fprintf('============================================================\n');


    ResultsPerGroup = struct();


    for j = 1:numel(sub_list)

        subj = sub_list{j};

        fprintf('\n');
        fprintf('>>> Processing %s for %s\n', subj, gname);


        %% ----------------------------------------------------------------
        % Construct subject-level analysis inputs
        % -----------------------------------------------------------------

        [~, ...
         ~, ...
         X1seed_train, ...
         X2seed_train, ...
         Xseed_test, ...
         Wg_0, ...
         TR, ...
         rseedtrain_vals, ...
         rseedtest_vals] = ...
            SubjectSeedfMRI( ...
                dataRoot, ...
                subj, ...
                vol_size_3d, ...
                atlas_file, ...
                atlas_data, ...
                seed_labels);


        %% ----------------------------------------------------------------
        % Store inputs
        % -----------------------------------------------------------------

        ResultsPerGroup(j).subject = subj;

        ResultsPerGroup(j).X1seed_train = X1seed_train;
        ResultsPerGroup(j).X2seed_train = X2seed_train;

        ResultsPerGroup(j).Wg_0 = Wg_0;

        ResultsPerGroup(j).Xseed_test = Xseed_test;

        ResultsPerGroup(j).Kcomp = Kcomp;
        ResultsPerGroup(j).TR = TR;

        ResultsPerGroup(j).seed_labels = seed_labels;

        ResultsPerGroup(j).rseedtrain_vals = rseedtrain_vals;
        ResultsPerGroup(j).rseedtest_vals = rseedtest_vals;


        fprintf('>>> %s / %s completed.\n', subj, gname);

    end


    % Store all subjects for the current ROI
    rsfMRIInputs.(gname) = ResultsPerGroup;

end


%% ========================================================================
% Save prepared inputs
% =========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Saving prepared rs-fMRI inputs...\n');
fprintf('============================================================\n');

save(outputFile, 'rsfMRIInputs', '-v7.3');


fprintf('\nPreparation completed successfully.\n');
fprintf('Output file:\n%s\n', outputFile);

fprintf('\nROIs included:\n');

for id_l = 1:numel(group_name)

    gname = group_name{id_l};

    fprintf('  %-5s : %d subjects\n', ...
        gname, numel(rsfMRIInputs.(gname)));

end

fprintf('\nDone.\n');



%% ========================================================================
% Local function: subject-level rs-fMRI input construction
% =========================================================================

function [resliced_filename, ...
          seed_indices, ...
          X1seed_train, ...
          X2seed_train, ...
          Xseed_test, ...
          Wg_0, ...
          TR, ...
          rseedtrain_vals, ...
          rseedtest_vals] = ...
          SubjectSeedfMRI( ...
              dataRoot, ...
              subj, ...
              vol_size_3d, ...
              atlas_file, ...
              atlas_data, ...
              seed_labels)


%% ========================================================================
% Subject paths
% =========================================================================

subj_dir = fullfile(dataRoot, subj);
func_dir = fullfile(subj_dir, 'func');

if ~exist(func_dir, 'dir')
    error('[%s] Functional-data directory not found:\n%s', ...
        subj, func_dir);
end


%% ========================================================================
% Locate the dataset-provided preprocessed functional image
% =========================================================================

preproc_filename = ...
    sprintf('%s_task-rest_space-MNI305_preproc.nii', subj);

preproc_fmri_file = ...
    fullfile(func_dir, preproc_filename);

if ~isfile(preproc_fmri_file)
    error('[%s] Preprocessed fMRI file not found:\n%s', ...
        subj, preproc_fmri_file);
end


%% ========================================================================
% Reslice fMRI to the Neuromorphometrics atlas grid
% =========================================================================
%
% The functional image provided with ds005747 is already preprocessed.
% Here it is resliced to the voxel grid of the Neuromorphometrics atlas
% so that the atlas-defined ROI voxel indices can be applied directly to
% the functional data.
%

resliced_filename = ...
    sprintf('r_%s_task-rest_space-MNI305_preproc.nii', subj);

resliced_fmri_file = ...
    fullfile(func_dir, resliced_filename);


if ~isfile(resliced_fmri_file)

    fprintf(['[%s] Reslicing preprocessed fMRI to the ', ...
             'Neuromorphometrics atlas grid...\n'], subj);

    matlabbatch = [];

    % Reference image: Neuromorphometrics atlas
    matlabbatch{1}.spm.spatial.coreg.write.ref = ...
        {atlas_file};

    % Source image: dataset-provided preprocessed fMRI
    matlabbatch{1}.spm.spatial.coreg.write.source = ...
        {preproc_fmri_file};

    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r_';

    spm_jobman('run', matlabbatch);


    if ~isfile(resliced_fmri_file)
        error('[%s] Resliced fMRI file was not created:\n%s', ...
            subj, resliced_fmri_file);
    end

else

    fprintf('[%s] Reusing existing atlas-grid resliced fMRI.\n', subj);

end


%% ========================================================================
% Load the atlas-grid resliced fMRI
% =========================================================================

fmri_vols = spm_vol(resliced_fmri_file);

T = numel(fmri_vols);
vol_size = prod(vol_size_3d);

% ===  Find seed voxel indices ===
seed_mask = ismember(atlas_data, seed_labels);
seed_indices = find(seed_mask); % seed voxel indices in MNI space
Label_timeseries = zeros(T, numel(seed_indices));

for t = 1:T
    vol_data = spm_read_vols(fmri_vols(t));
    Label_timeseries(t, :) = vol_data(seed_indices);
end

fprintf('=====DMN_timeseries completed====\n');

% Split run1 into training and testing in time
T1 = floor(T / 3);
TR = (10 * 60) / T1;
T_half = T1 / 2;

Label_train = Label_timeseries(1:T_half, :);     %[128 x V]
Label_test  = Label_timeseries(T_half+1:T1, :);  %[128 x V]


% --- Build 4D volume for (training) and testing with intrested seed voxel only---
fMRI_train_flat = zeros(vol_size, T_half);
fMRI_train_flat(seed_indices, :) = Label_train';
fMRI_train_4D = reshape(fMRI_train_flat, [vol_size_3d T_half]);


X1seed_train = Label_train'; %[V_seed x 128]

% === Compute seed-based r-values for training ===
% Z-score each voxel time series (row-wise z-score)
voxel_z = zscore(X1seed_train, 0, 2);          % [V_seed x T]
% Average seed time series across voxels, then z-score
seedmean_z = zscore(mean(voxel_z, 1));         % [1 x T]
% Correlation (dot product of z-scores)
rseedtrain_vals = (voxel_z * seedmean_z') / (T_half - 1);  % [V_seed x 1]

V=size(X1seed_train, 1);
Z0 = [ones(V,1), rseedtrain_vals];   % Design matrix
wg0_vec = zeros(T_half, 1);

for t = 1:T_half
    x_t = X1seed_train(:, t);   % voxel values at time t
    beta_t = (Z0' * Z0) \ (Z0' * x_t);
    wg0_vec(t) = norm(beta_t);
end

wg0_vec = T * wg0_vec / sum(wg0_vec);
Wg_0 = diag(wg0_vec);

[~, Z0_idx]=sort(rseedtrain_vals, 'descend');
Z0_label_idx=seed_indices(Z0_idx(1:floor(V*0.2)));
ZO_mask=zeros(vol_size,1);
ZO_mask(Z0_label_idx,:)=1;
ZO_3D_mask=reshape(ZO_mask,vol_size_3d);



% === Compute the second-view for training ===

[X2_4D, ~] = BuildActiveCCAViews(fMRI_train_4D, seed_mask, ZO_3D_mask);
X2=reshape(X2_4D,vol_size,T_half);
X2seed_train=X2(seed_indices,:);
Xseed_test=Label_test';

% === Compute seed-based r-values for test ===
voxeltest_z = zscore(Xseed_test, 0, 2);         % [V x T]
seedtest_mean_z = zscore(mean(voxeltest_z, 1)); % [1 x T]
rseedtest_vals = (voxeltest_z * seedtest_mean_z') / (T_half - 1);  % [V x 1]


end


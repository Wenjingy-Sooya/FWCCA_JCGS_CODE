% =========================================================================
% RunGLMwithThreeConditions.m
%
% Subject-level SPM first-level GLM analysis for the OpenNeuro Emotional
% Music dataset (ds000171).
%
% This script fits a first-level GLM to the preprocessed trainrun fMRI data
% for one subject using three experimental conditions:
%
%    1. tones
%    2. positive_music
%    3. negative_music
%
% INPUTS
% ------
% 1. Preprocessed training fMRI:
%
%    Data/<group>/<subject>/func/
%        swar<subject>_task-music_trainrun_bold.nii
%
% 2. Subject-specific event information for the three experimental
%    conditions:
%
%        tones
%        positive_music
%        negative_music
%
% The GLM is fitted using the trainrun only. The independent testrun is
% not used for first-level model estimation.
%
% MODEL
% -----
% The first-level SPM model contains separate regressors for the three
% experimental conditions. After model specification and estimation,
% condition-specific contrasts are computed for the subsequent FWCCA
% analysis.
%
% OUTPUTS
% -------
% Subject-specific first-level results are saved under:
%
%    Data/<group>/<subject>/first_level_analysis/
%
% including:
%
%    SPM.mat
%    mask.nii
%    spmT_*.nii
%
% The resulting first-level GLM outputs are subsequently used to construct
% the positive and negative second views in
% 03_SecondViewSensitivityAnalysis.
%
% IMPORTANT
% ---------
% The first-level GLM is estimated from the trainrun only. In particular,
% the GLM-derived statistics used for second-view construction are obtained
% without using the held-out testrun.
% =========================================================================

clear;
close all;
clc;
clear; close all; clc;

%% ============================================================
% Project and SPM setup
% ============================================================

glmRoot     = fileparts(mfilename('fullpath'));
datasetRoot = fileparts(glmRoot);
projectRoot = fileparts(datasetRoot);

spm_dir = fullfile(projectRoot, 'spm_25');

if ~exist(spm_dir, 'dir')
    error('SPM directory does not exist:\n%s', spm_dir);
end

addpath(spm_dir);

spm('defaults', 'fmri');
spm_jobman('initcfg');


%% ============================================================
% Data paths and analysis settings
% ============================================================

group = 'MDD';   % Options: 'ND' or 'MDD'

switch group

    case 'ND'

        sub_list = { ...
            'sub-control01', ...
            'sub-control04', ...
            'sub-control06', ...
            'sub-control07', ...
            'sub-control09', ...
            'sub-control11', ...
            'sub-control15', ...
            'sub-control16', ...
            'sub-control17', ...
            'sub-control18'};

    case 'MDD'

        sub_list = { ...
            'sub-mdd02', ...
            'sub-mdd03', ...
            'sub-mdd05', ...
            'sub-mdd07', ...
            'sub-mdd09', ...
            'sub-mdd10', ...
            'sub-mdd14', ...
            'sub-mdd15', ...
            'sub-mdd18'};

    otherwise

        error('Unknown group: %s', group);

end

dataRoot = fullfile(datasetRoot, 'Data', group);

% The first-level GLM is fitted using the training run.
trainrun = 'trainrun';


%% ============================================================
% First-level GLM
% ============================================================

for s = 1:numel(sub_list)

    subj = sub_list{s};

    fprintf('\n============================================================\n');
    fprintf('Running first-level GLM for %s\n', subj);
    fprintf('============================================================\n');

    %% --------------------------------------------------------
    % Subject directories
    % ---------------------------------------------------------

    subj_dir = fullfile(dataRoot, subj);
    func_dir = fullfile(subj_dir, 'func');

    % Use a separate directory for the checking run so that the
    % previously saved first-level analysis is not overwritten.
    outdir = fullfile(subj_dir, 'first_level_analysis_checking');

    if ~exist(subj_dir, 'dir')
        warning('%s: subject directory does not exist. Skipping.', subj);
        continue;
    end

    if ~exist(func_dir, 'dir')
        warning('%s: functional directory does not exist. Skipping.', subj);
        continue;
    end

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end


    %% --------------------------------------------------------
    % Locate preprocessed functional data
    % ---------------------------------------------------------

    func_pattern = ...
        ['^swar' subj '_task-music_' trainrun '_bold\.nii$'];

    nii_file = spm_select( ...
        'FPList', ...
        func_dir, ...
        func_pattern);

    motion_pattern = ...
        ['^rp_' subj '_task-music_' trainrun '_bold\.txt$'];

    motion_file = spm_select( ...
        'FPList', ...
        func_dir, ...
        motion_pattern);


    %% --------------------------------------------------------
    % Check required inputs
    % ---------------------------------------------------------

    if isempty(nii_file)

        warning(['%s: preprocessed training-run functional image ' ...
                 'was not found. Skipping.'], subj);

        continue;

    end

    if isempty(motion_file)

        warning(['%s: motion parameter file was not found. ' ...
                 'Skipping.'], subj);

        continue;

    end


    %% --------------------------------------------------------
    % Expand the 4-D NIfTI into individual scans for SPM
    % ---------------------------------------------------------

    scans = cellstr( ...
        spm_select( ...
        'ExtFPList', ...
        func_dir, ...
        func_pattern, ...
        Inf));

    if isempty(scans)

        warning('%s: no functional volumes were found. Skipping.', subj);
        continue;

    end

    fprintf('Functional image:\n%s\n', nii_file);
    fprintf('Motion parameters:\n%s\n', motion_file);
    fprintf('Number of scans: %d\n', numel(scans));
    fprintf('Output directory:\n%s\n\n', outdir);


    %% ========================================================
    % STEP 1: GLM specification
    % =========================================================

    matlabbatch = {};

    matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir};

    % Timing
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 3;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    % Functional scans
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;


    %% --------------------------------------------------------
    % Experimental conditions
    % ---------------------------------------------------------

    cond_names = { ...
        'tones', ...
        'positive_music', ...
        'negative_music'};

    cond_onsets = { ...
        [0, 70.5, 139.5, 208.5, 277.5], ...  % tones
        [105, 243], ...                       % positive_music
        [36, 174] ...                         % negative_music
        };

    cond_durations = { ...
        repmat(31.5, 1, 5), ...   % tones
        repmat(31.5, 1, 2), ...   % positive_music
        repmat(31.5, 1, 2) ...    % negative_music
        };

    for c = 1:numel(cond_names)

        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).name = ...
            cond_names{c};

        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).onset = ...
            cond_onsets{c};

        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).duration = ...
            cond_durations{c};

        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).tmod = 0;

        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).pmod = ...
            struct('name', {}, 'param', {}, 'poly', {});

        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).orth = 1;

    end


    %% --------------------------------------------------------
    % Motion regressors and model settings
    % ---------------------------------------------------------

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = ...
        struct('name', {}, 'val', {});

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = ...
        {motion_file};

    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

    matlabbatch{1}.spm.stats.fmri_spec.fact = ...
        struct('name', {}, 'levels', {});

    % Canonical HRF without derivatives
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];

    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};

    % Temporal autocorrelation model
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';


    %% ========================================================
    % STEP 2: Estimate GLM
    % =========================================================

    matlabbatch{2}.spm.stats.fmri_est.spmmat = ...
        {fullfile(outdir, 'SPM.mat')};

    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;

    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


    %% ========================================================
    % STEP 3: Define contrasts
    % =========================================================

    matlabbatch{3}.spm.stats.con.spmmat = ...
        {fullfile(outdir, 'SPM.mat')};


    % Contrast 1: tones
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ...
        'tones';

    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = ...
        [1 0 0];

    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = ...
        'none';


    % Contrast 2: positive music
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ...
        'positive_music';

    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = ...
        [0 1 0];

    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = ...
        'none';


    % Contrast 3: negative music
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = ...
        'negative_music';

    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = ...
        [0 0 1];

    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = ...
        'none';

    matlabbatch{3}.spm.stats.con.delete = 0;


    %% ========================================================
    % Save batch and run SPM
    % =========================================================

    batch_file = fullfile( ...
        outdir, ...
        'firstlevel_GLM_batch.mat');

    save(batch_file, 'matlabbatch');

    spm_jobman('run', matlabbatch);


    %% --------------------------------------------------------
    % Check expected outputs
    % ---------------------------------------------------------

    expected_outputs = { ...
        'SPM.mat', ...
        'spmT_0001.nii', ...
        'spmT_0002.nii', ...
        'spmT_0003.nii'};

    all_outputs_exist = true;

    for j = 1:numel(expected_outputs)

        current_file = fullfile(outdir, expected_outputs{j});

        if ~isfile(current_file)

            warning('%s: expected output not found: %s', ...
                subj, expected_outputs{j});

            all_outputs_exist = false;

        end

    end

    if all_outputs_exist

        fprintf('\n--- %s completed successfully ---\n', subj);
        fprintf('spmT_0001.nii : tones\n');
        fprintf('spmT_0002.nii : positive_music\n');
        fprintf('spmT_0003.nii : negative_music\n');

    else

        fprintf('\n--- %s completed, but some outputs are missing ---\n', ...
            subj);

    end

end

fprintf('\n============================================================\n');
fprintf('First-level GLM analysis finished for group: %s\n', group);
fprintf('============================================================\n');
% =========================================================================
% TwophasePreprocessing.m
%
% Two-phase SPM preprocessing pipeline for the OpenNeuro Emotional Music
% dataset (ds000171).
%
% This script preprocesses the anatomical and functional MRI data for each
% subject. The trainrun and testrun functional images are processed
% separately so that the two runs remain independent for the subsequent
% training and held-out testing analyses.
%
% PREPROCESSING WORKFLOW
% ----------------------
% Phase 1:
%
%    Realign
%       ->
%    Slice Timing
%       ->
%    Coregister
%       ->
%    Segment
%
% Phase 2:
%
%    Normalize
%       ->
%    Smooth
%
% INPUTS
% ------
% Subject-specific MRI data are read from:
%
%    Data/<group>/<subject>/
%
% with the following structure:
%
%    anat/
%        <subject>_T1w.nii
%
%    func/
%        <subject>_task-music_trainrun_bold.nii
%        <subject>_task-music_testrun_bold.nii
%
% where <group> is either MDD or ND.
%
% PROCESSING
% ----------
% The preprocessing procedure is applied separately to the trainrun and
% testrun functional images.
%
% Phase 1 performs:
%
%    1. Functional-image realignment
%    2. Slice-timing correction
%    3. Coregistration of functional and anatomical images
%    4. Segmentation of the anatomical image
%
% Phase 2 performs:
%
%    5. Spatial normalization
%    6. Spatial smoothing
%
% OUTPUTS
% -------
% SPM-generated preprocessing outputs are stored in the corresponding
% subject-specific anat/ and func/ directories.
%
% The final smoothed, normalized, slice-time-corrected, and realigned
% functional images follow the naming convention:
%
%    swar<subject>_task-music_trainrun_bold.nii
%    swar<subject>_task-music_testrun_bold.nii
%
% These preprocessed functional images are subsequently used in the
% first-level GLM and FWCCA analyses.
%
% IMPORTANT
% ---------
% The trainrun and testrun are preprocessed separately. The trainrun is
% subsequently used for model construction and hyperparameter selection,
% whereas the testrun is reserved for independent held-out evaluation.
% =========================================================================

clear;
close all;
clc;


%% ============================================================
% Project and SPM setup
% ============================================================

preprocessingRoot = fileparts(mfilename('fullpath'));
datasetRoot = fileparts(preprocessingRoot);
projectRoot = fileparts(datasetRoot);


spm_dir = fullfile(projectRoot, 'spm_25');

if ~exist(spm_dir, 'dir')
    error('SPM directory does not exist:\n%s', spm_dir);
end

addpath(spm_dir);
spm('defaults', 'fmri');
spm_jobman('initcfg');

%% ============================================================
% Data paths and preprocessing settings
% ============================================================

group = 'MDD';   % Options: 'ND' or 'MDD'

switch group
    case 'ND'
        sub_list = { ...
            'sub-control01', 'sub-control04', 'sub-control06', ...
            'sub-control07', 'sub-control09', 'sub-control11', ...
            'sub-control15', 'sub-control16', 'sub-control17', ...
            'sub-control18'};

    case 'MDD'
        %{
        sub_list = { ...
            'sub-mdd02', 'sub-mdd03', 'sub-mdd05', ...
            'sub-mdd07', 'sub-mdd09', 'sub-mdd10', ...
            'sub-mdd14', 'sub-mdd15', 'sub-mdd18'};
        %}
         sub_list={'sub-mdd02'};

    otherwise
        error('Unknown group: %s', group);
end

dataRoot = fullfile(datasetRoot, 'Data', group);
func_runs = {'trainrun','testrun'};     

for i = 1:numel(sub_list)

    subj = sub_list{i};
    subj_dir = fullfile(dataRoot, subj);
    anat_dir = fullfile(subj_dir, 'anat');
    func_dir = fullfile(subj_dir, 'func');

    % Get anatomical T1 path
    anat = spm_select('FPList', anat_dir, ['^' subj '_T1w\.nii$']);

    for r = 1:numel(func_runs)
        runname = func_runs{r};
        run_prefix = [subj '_task-music_' runname];
        func_path = fullfile(func_dir, [run_prefix '_bold.nii']);

        if ~isfile(func_path)
            warning('Missing functional file: %s — skipping.', func_path);
            continue;
        end

        %% === Phase 1 ===
        fprintf('\n[Phase 1] Starting %s — %s\n', subj, runname);
        matlabbatch = {};

        % Step 1: Realign
        matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = {func_path};
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.95;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 1.5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

        % Step 2: Slice Timing
        matlabbatch{2}.spm.temporal.st.scans{1} = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', ...
            substruct('.', 'val', '{}',{1}, '.', 'val', '{}',{1}, '.', 'val', '{}',{1}, '.', 'val', '{}',{1}), ...
            substruct('.', 'sess', '()', {1}, '.', 'cfiles'));
        matlabbatch{2}.spm.temporal.st.nslices = 50;
        matlabbatch{2}.spm.temporal.st.tr = 3.0;
        matlabbatch{2}.spm.temporal.st.ta = 3.0 - (3.0 / 50);
        matlabbatch{2}.spm.temporal.st.so = [1:2:50, 2:2:50];
        matlabbatch{2}.spm.temporal.st.refslice = 25;
        matlabbatch{2}.spm.temporal.st.prefix = 'ar';

        % Step 3: Coregister
        matlabbatch{3}.spm.spatial.coreg.estwrite.ref = cellstr(anat);
        matlabbatch{3}.spm.spatial.coreg.estwrite.source(1) = cfg_dep( ...
            'Realign: Estimate & Reslice: Mean Image', ...
            substruct('.', 'val', '{}',{1}, '.', 'val', '{}',{1}, '.', 'val', '{}',{1}, '.', 'val', '{}',{1}), ...
            substruct('.', 'rmean'));
        matlabbatch{3}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.prefix = 'c';

        % Step 4: Segment
        matlabbatch{4}.spm.spatial.preproc.channel.vols = cellstr(anat);
        matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.0001;
        matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1];
        ngausorder = [1 1 2 3 4 2];
        for k = 1:6
            matlabbatch{4}.spm.spatial.preproc.tissue(k).tpm = {[fullfile(spm('Dir'),'tpm','TPM.nii') ',' num2str(k)]};
            matlabbatch{4}.spm.spatial.preproc.tissue(k).ngaus = ngausorder(k);
            matlabbatch{4}.spm.spatial.preproc.tissue(k).native = [1 0];
            matlabbatch{4}.spm.spatial.preproc.tissue(k).warped = [0 0];
        end
        matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0 0.1 0.01 0.04];
        matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];

        save(fullfile(func_dir, ['phase1_' runname '_batch' '.mat']), 'matlabbatch');
        spm_jobman('run', matlabbatch);

        %% === Phase 2 ===
        fprintf('\n[Phase 2] Starting normalization & smoothing for %s — %s\n', subj, runname);
        matlabbatch = {};

        def_field = spm_select('FPList', fullfile(subj_dir, 'anat'), '^y_.*\.nii$');
        if isempty(def_field)
            warning('No deformation field found for %s — skipping Phase 2.', subj);
            continue;
        end

        ar_func = spm_select('FPList', func_dir, ['^ar' run_prefix '.*\.nii$']);

        if isempty(ar_func)
            warning('No slice-timed images found for %s — skipping Phase 2.', run_prefix);
            continue;
        end

        % Step 5: Normalize
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(def_field);
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(ar_func);
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        % Step 6: Smooth
        matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep( ...
            'Normalise: Write: Normalised Images (Subj 1)', ...
            substruct('.', 'val', '{}',{1}, '.', 'val', '{}',{1}, '.', 'val', '{}',{1}, '.', 'val', '{}',{1}), ...
            substruct('()', {1}, '.', 'files'));
        matlabbatch{2}.spm.spatial.smooth.fwhm = [4 4 4];
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's';

        save(fullfile(func_dir, ['phase2_' runname '_batch' '.mat']), 'matlabbatch');
        spm_jobman('run', matlabbatch);
        
    end
end
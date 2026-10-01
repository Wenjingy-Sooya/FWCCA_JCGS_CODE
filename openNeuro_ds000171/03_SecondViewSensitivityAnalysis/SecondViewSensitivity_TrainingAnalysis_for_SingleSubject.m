% =========================================================================
% SecondViewSensitivity_TrainingAnalysis_for_SingleSubject.m
%
% Training analysis for the FWCCA second-view sensitivity study using the
% OpenNeuro Emotional Music dataset (ds000171).
%
% This script uses the trainrun data from one subject to construct positive
% and negative second views at six percentile thresholds and to fit CCA and
% G-FWCCA for each second-view definition.
%
% The script is provided as a single-subject example using sub-mdd02 from
% the MDD group. The same procedure can be applied to other MDD or ND
% subjects by changing groupName and subj below.
%
% INPUTS
% ------
% 1. Preprocessed training fMRI:
%    Data/<group>/<subject>/func/
%        swar<subject>_task-music_trainrun_bold.nii
%
% 2. Subject-specific first-level SPM model:
%    Data/<group>/<subject>/first_level_analysis/SPM.mat
%
% 3. First-level GLM t-maps:
%    Data/<group>/<subject>/first_level_analysis/
%        spmT_0002.nii    positive-music contrast
%        spmT_0003.nii    negative-music contrast
%
% 4. First-level analysis mask:
%    Data/<group>/<subject>/first_level_analysis/mask.nii
%
% 5. Subject-specific tissue probability maps:
%    Data/<group>/<subject>/anat/
%
%    These are resliced to the first-level analysis mask and used to
%    construct the cleaned anatomical analysis mask.
%
% SECOND-VIEW DEFINITIONS
% -----------------------
% Positive and negative second views are constructed separately from the
% corresponding training-derived GLM t-maps.
%
% Six percentile thresholds are considered:
%
%    p90 -> top 10%
%    p85 -> top 15%
%    p80 -> top 20%
%    p75 -> top 25%
%    p70 -> top 30%
%    p65 -> top 35%
%
% For each threshold, the selected voxels define a candidate mask. The
% candidate mask is then used by BuildActiveCCAViews to construct the
% corresponding second view, including the spatial dilation implemented
% by that function.
%
% GLOBAL FEATURE-WEIGHT INITIALIZATION
% ------------------------------------
% The initial global feature weights Wg_0 are generated using a fixed
% random seed and are intentionally independent of the GLM t-maps. The
% same initialization is used for all second-view thresholds and for both
% positive and negative second-view analyses.
%
% METHODS FITTED IN THIS SCRIPT
% -----------------------------
% For each second-view definition and threshold:
%
%    1. CCA
%    2. G-FWCCA
%
% CCA uses identity global weights and identity local kernels.
% G-FWCCA uses Wg_0 as the global feature-weight matrix and identity local
% kernels.
%
% LG-FWCCA and GL-FWCCA are not fitted in this script. Their local-kernel
% hyperparameters are selected subsequently using
% SecondViewSensitivity_LG_GLFWCCA_CV_for_SingleSubject.m.
%
% OUTPUTS
% -------
% Results are saved under:
%
%    03_SecondViewSensitivityAnalysis/
%        SubjectResults/<group>/<subject>/
%
% Two MAT files are generated:
%
% 1. <subject>_SecondViewInputs.mat
%
%    Contains AllTrainInputs, including the training analysis mask, X1,
%    Wg_0, analysis settings, and the threshold-specific second views X2.
%    This file is subsequently used for LG-FWCCA and GL-FWCCA
%    cross-validation and for held-out testing.
%
% 2. <subject>_CCA_GFWCCA_TrainingResults.mat
%
%    Contains CCA_GFWCCA_TrainingResults with the fitted CCA and G-FWCCA
%    results for the positive and negative second views at all six
%    thresholds.
%
% IMPORTANT
% ---------
% All second-view construction and model fitting in this script use the
% trainrun only. The independent testrun is not used at this stage.
%
% Example subject:
%    Group   : MDD
%    Subject : sub-mdd02
%    Run     : trainrun
% =========================================================================

clear; close all; clc;

%% ============================================================
% Project paths
% ============================================================
analysisRoot = fileparts(mfilename('fullpath'));
datasetRoot  = fileparts(analysisRoot);
projectRoot  = fileparts(datasetRoot);

spm_dir       = fullfile(projectRoot, 'spm_25');
functions_dir = fullfile(datasetRoot, 'Functions');

if ~exist(spm_dir, 'dir')
    error('SPM directory does not exist:\n%s', spm_dir);
end

if ~exist(functions_dir, 'dir')
    error('Functions directory does not exist:\n%s', functions_dir);
end

%% ============================================================
% SPM and function setup
% ============================================================
addpath(spm_dir);
addpath(functions_dir);

spm('defaults', 'fmri');
spm_jobman('initcfg');

atlas = spm_atlas('load', 'Neuromorphometrics');

%% ============================================================
% Subject settings
%
% sub-mdd02 from the MDD group is used as an example.
% ============================================================

groupName = 'MDD';
subj = 'sub-mdd02';
currentrun = 'trainrun';

dataRoot = fullfile(datasetRoot, 'Data', groupName);

%% ============================================================
% Analysis settings
% ============================================================

viewTypeList = {'positive', 'negative'};
Kcomp_max = 5;
svdTol = 1e-2;
ridgeTol = 1e-5;
percentile_list = [90 85 80 75 70 65];
nThresholds = numel(percentile_list);

fprintf('\n');
fprintf('============================================================\n');
fprintf('Processing subject %s\n',subj);
fprintf('============================================================\n');

%% --------------------------------------------------------
% Subject paths
% --------------------------------------------------------
subj_dir = fullfile(dataRoot,subj);
func_dir = fullfile(subj_dir,'func');
anat_dir = fullfile(subj_dir,'anat');
spm_path = fullfile(subj_dir,'first_level_analysis','SPM.mat');

if ~exist(subj_dir,'dir')
    error('Subject directory not found:\n%s',subj_dir);
end

if ~exist(spm_path,'file')
    error('SPM.mat not found:\n%s',spm_path);
end

%% --------------------------------------------------------
% Load functional data
% --------------------------------------------------------
functional_pattern = ['^swar' subj '_task-music_' currentrun '_bold\.nii$'];
files = spm_select('FPList',func_dir,functional_pattern);

if isempty(files)
    error('Cannot find preprocessed fMRI for %s.',subj);
end

nii_path = strtrim(files(1,:));
fprintf('Functional image:\n%s\n',nii_path);

%% --------------------------------------------------------
% Load SPM model
% --------------------------------------------------------
load(spm_path,'SPM');
analysis_dir = fileparts(spm_path);
fprintf('Stored SPM.swd:\n%s\n',SPM.swd);
fprintf('Using current analysis directory:\n%s\n',analysis_dir);
TR = SPM.xY.RT;
fprintf('TR = %.3f sec\n',TR);

%% --------------------------------------------------------
% Drift removal and temporal prewhitening
% --------------------------------------------------------
Y = get_driftremoved_fMRI(spm_path,nii_path);
V = full(SPM.xVi.V);
if size(V,1)~=size(Y,2)
    error('Temporal covariance dimension mismatch.');
end

Vw = sqrtm(inv(V));
Y = Y*Vw;
if ~isreal(Y)
    warning('Complex values produced during whitening. Taking real part.');
    Y = real(Y);
end

T = size(Y,2);
fprintf('Time points = %d\n',T);

%% --------------------------------------------------------
% Construct anatomical analysis mask
% --------------------------------------------------------
mask_file = fullfile(analysis_dir,'mask.nii');

if ~exist(mask_file,'file')
    error('Cannot find analysis mask:\n%s',mask_file);
end

rc1_file = reslice_to_mask(mask_file,anat_dir,subj,'c1');
rc2_file = reslice_to_mask(mask_file,anat_dir,subj,'c2');
rc3_file = reslice_to_mask(mask_file,anat_dir,subj,'c3');

opts = struct();
opts.exclude_labels = {...
    'Brain Stem',...
    '4th Ventricle',...
    '3rd Ventricle',...
    'Left Lateral Ventricle',...
    'Right Lateral Ventricle'};

[final_mask,idx_mask,idx_non,info] = ...
    generate_clean_mask(mask_file,...
                        rc1_file,...
                        rc2_file,...
                        rc3_file,...
                        opts);
vol_size_3d = size(final_mask);
num_vox = prod(vol_size_3d);
nMaskVoxels = numel(idx_mask);
fprintf('Analysis mask voxels = %d\n',nMaskVoxels);

%% --------------------------------------------------------
% Rebuild 4D image
% --------------------------------------------------------
Y_4D = nan([vol_size_3d T]);

for t = 1:T
    tmp = zeros(num_vox,1);
    tmp(idx_mask) = Y(idx_mask,t);
    Y_4D(:,:,:,t) = reshape(tmp,vol_size_3d);
end

%% --------------------------------------------------------
% Load GLM t-maps
% --------------------------------------------------------
tmap_pos_file = fullfile(analysis_dir,'spmT_0002.nii');
tmap_neg_file = fullfile(analysis_dir,'spmT_0003.nii');

if ~exist(tmap_pos_file,'file')
    error('Cannot find positive t-map:\n%s',tmap_pos_file);
end

if ~exist(tmap_neg_file,'file')
    error('Cannot find negative t-map:\n%s',tmap_neg_file);
end

tmap_pos = spm_read_vols(spm_vol(tmap_pos_file));
tmap_neg = spm_read_vols(spm_vol(tmap_neg_file));
tmap_pos_masked = tmap_pos(idx_mask);
tmap_neg_masked = tmap_neg(idx_mask);
validPositive = isfinite(tmap_pos_masked);
validNegative = isfinite(tmap_neg_masked);

%% --------------------------------------------------------
% Random fixed initialization of global feature weights
%
% W0 is intentionally independent of the GLM t-maps.
% The same random initialization is used for every threshold.
%% --------------------------------------------------------
rng(0);
wg0_vec = rand(T,1);
wg0_vec = T*wg0_vec/sum(wg0_vec);
Wg_0 = diag(wg0_vec);
X1Mask = Y(idx_mask,:);

%% ============================================================
% Run second-view sensitivity analysis
%% ============================================================
nViews = numel(viewTypeList);

AllTrainInputs  = cell(nViews, 1);
CCA_GFWCCA_TrainingResults = cell(nViews, 1);

for viewIdx = 1:nViews

     viewType = viewTypeList{viewIdx};
     fprintf('\n');
     fprintf('============================================================\n');
     fprintf('Running %s second-view analysis\n', viewType);
     fprintf('============================================================\n');

    %% --------------------------------------------------------
    % Select contrast map
    %% --------------------------------------------------------

    switch lower(viewType)

        case 'positive'

            contrastMap = tmap_pos_masked;
            validMask   = validPositive;

        case 'negative'

            contrastMap = tmap_neg_masked;
            validMask   = validNegative;

        otherwise

            error('Unknown viewType: %s', viewType);

    end

    %% --------------------------------------------------------
    % Initialize structures
    %% --------------------------------------------------------
    trainInput               = struct();
    trainOutput              = struct();
    trainInput.viewType      = viewType;
    trainOutput.viewType     = viewType;
    trainInput.subject       = subj;
    trainInput.currentRun    = currentrun;
    trainInput.X1Mask        = X1Mask;
    trainInput.idx_mask      = idx_mask;
    trainInput.idx_non       = idx_non;
    trainInput.final_mask    = final_mask;
    trainInput.maskInfo      = info;
    trainInput.volumeSize    = vol_size_3d;
    trainInput.nMaskVoxels   = nMaskVoxels;
    trainInput.TR            = TR;
    trainInput.T             = T;
    trainInput.Kcomp_max     = Kcomp_max;
    trainInput.svdTol        = svdTol;
    trainInput.Wg_0          = Wg_0;
    trainInput.wg0_vec       = wg0_vec;
    %% --------------------------------------------------------
    % Threshold template
    %% --------------------------------------------------------
    thresholdTemplate = struct( ...
        'percentile', [], ...
        'topProportion', [], ...
        'thresholdValue', [], ...
        'nSelectedBeforeDilation', [], ...
        'proportionBeforeDilation', [], ...
        'selectedMaskIndicesWithinAnalysisMask', [], ...
        'selectedMaskLinearIndicesFullVolume', [], ...
        'selectedMaskBeforeDilation', [], ...
        'nSelectedAfterDilation', [], ...
        'proportionAfterDilation', [], ...
        'activeMaskDilated', [], ...
        'X2Mask', []);

    trainInput.secondViewThresholds = repmat(thresholdTemplate, nThresholds, 1);
    %% --------------------------------------------------------
    % Output template
    %% --------------------------------------------------------
    resultTemplate = struct( ...
        'percentile', [], ...
        'topProportion', [], ...
        'resCCA_T', [], ...
        'resG_FWCCA_T', []);

    trainOutput.subject    = subj;
    trainOutput.currentRun = currentrun;
    trainOutput.viewType   = viewType;
    trainOutput.X1Mask     = X1Mask;
    trainOutput.final_mask = final_mask;
    trainOutput.idx_mask   = idx_mask;
    trainOutput.results    = repmat(resultTemplate, nThresholds, 1);

    %% ========================================================
    % Loop over threshold levels
    %% ========================================================

    for id_p = 1:nThresholds

        percentileValue = percentile_list(id_p);
        topProportion   = 100 - percentileValue;

        fprintf('\n');
        fprintf('%s: constructing %s second view using top %d%% (percentile = %d)\n', ...
            subj, viewType, topProportion, percentileValue);

        %% ----------------------------------------------------
        % Candidate region
        %% ----------------------------------------------------

        thresholdValue = prctile( ...
            contrastMap(validMask), ...
            percentileValue);

        candidate_idx = find( ...
            isfinite(contrastMap) & ...
            contrastMap >= thresholdValue);

        nSelectedBeforeDilation = numel(candidate_idx);

        proportionBeforeDilation = ...
            100 * nSelectedBeforeDilation / nMaskVoxels;

        %% ----------------------------------------------------
        % Candidate mask
        %% ----------------------------------------------------
        selectedFullVolumeIndices = idx_mask(candidate_idx);
        candidateMask1D = zeros(num_vox,1);
        candidateMask1D(selectedFullVolumeIndices) = 1;
        candidateMask3D = reshape(candidateMask1D, vol_size_3d);
        %% ----------------------------------------------------
        % Construct second view
        %% ----------------------------------------------------

        [Y2_4D, activeMaskDilated] = BuildActiveCCAViews( ...
                Y_4D, ...
                final_mask, ...
                candidateMask3D);

        X2 = reshape(Y2_4D, num_vox, T);
        X2Mask = X2(idx_mask,:);

        if ~isequal(size(X2Mask), size(X1Mask))

              error('Second-view dimensions do not match X1Mask.');
        end

        %% ----------------------------------------------------
        % Summarize dilated mask
        %% ----------------------------------------------------

        if isempty(activeMaskDilated)

            nSelectedAfterDilation = NaN;
            proportionAfterDilation = NaN;

        else

            activeMaskDilated = logical(activeMaskDilated);

            nSelectedAfterDilation = ...
                sum(activeMaskDilated(idx_mask));

            proportionAfterDilation = ...
                100 * nSelectedAfterDilation / nMaskVoxels;

        end

        %% ----------------------------------------------------
        % Store threshold information
        %% ----------------------------------------------------

        trainInput.secondViewThresholds(id_p).percentile = percentileValue;
        trainInput.secondViewThresholds(id_p).topProportion = topProportion;
        trainInput.secondViewThresholds(id_p).thresholdValue = thresholdValue;
        trainInput.secondViewThresholds(id_p).nSelectedBeforeDilation = nSelectedBeforeDilation;
        trainInput.secondViewThresholds(id_p).proportionBeforeDilation = proportionBeforeDilation;
        trainInput.secondViewThresholds(id_p).selectedMaskIndicesWithinAnalysisMask = candidate_idx;
        trainInput.secondViewThresholds(id_p).selectedMaskLinearIndicesFullVolume = selectedFullVolumeIndices;
        trainInput.secondViewThresholds(id_p).selectedMaskBeforeDilation = candidateMask3D;
        trainInput.secondViewThresholds(id_p).nSelectedAfterDilation = nSelectedAfterDilation;
        trainInput.secondViewThresholds(id_p).proportionAfterDilation = proportionAfterDilation;
        trainInput.secondViewThresholds(id_p).activeMaskDilated = activeMaskDilated;
        trainInput.secondViewThresholds(id_p).X2Mask = X2Mask;

        
        fprintf('Before dilation: %d voxels (%.2f%% of analysis mask)\n',nSelectedBeforeDilation, proportionBeforeDilation);

        if isfinite(nSelectedAfterDilation)
            fprintf('After dilation: %d voxels (%.2f%% of analysis mask)\n', nSelectedAfterDilation,  proportionAfterDilation);
        end

        %% --------------------------------------------------------
        % Run CCA variants
        %% --------------------------------------------------------
        fprintf('Running CCA variants...\n');

        resCCA_T = Run_FWCCA_Family( ...
            X1Mask, X2Mask, ...
            eye(T), eye(T), ...
            eye(T), eye(T), ...
            Kcomp_max, 100, 1e-5, ...
            'prescale', ridgeTol, svdTol);

        resG_FWCCA_T = Run_FWCCA_Family( ...
            X1Mask, X2Mask, ...
            Wg_0, Wg_0, ...
            eye(T), eye(T), ...
            Kcomp_max, 100, 1e-5, ...
            'prescale', ridgeTol, svdTol);


        %% --------------------------------------------------------
        % Keep only required CCA / G-FWCCA fields
        %% --------------------------------------------------------

        fieldsToKeep = { 'Ak','Bk', ...
                         'Ak_raw', 'Bk_raw', ...
                         'rhoVec', 'X1_localMean', ...
                         'X2_localMean',  'mixOrder', ...
                         'effectiveWeight'};

        fieldsToRemove = setdiff( fieldnames(resCCA_T), ...
                                 fieldsToKeep);

        if ~isempty(fieldsToRemove)

          resCCA_T = rmfield( resCCA_T, fieldsToRemove);

        end

        fieldsToRemove = setdiff(fieldnames(resG_FWCCA_T), fieldsToKeep);

        if ~isempty(fieldsToRemove)

            resG_FWCCA_T = rmfield(  resG_FWCCA_T,  fieldsToRemove);

        end

     
        %% --------------------------------------------------------
        % Store results for this threshold
        %% --------------------------------------------------------
        trainOutput.results(id_p).percentile = percentileValue;
        trainOutput.results(id_p).topProportion = topProportion;
        trainOutput.results(id_p).resCCA_T = resCCA_T;
        trainOutput.results(id_p).resG_FWCCA_T = resG_FWCCA_T;

    end
      %% --------------------------------------------------------
      % Store this view
      %% --------------------------------------------------------

       AllTrainInputs{viewIdx}  = trainInput;
       CCA_GFWCCA_TrainingResults{viewIdx} = trainOutput;

     
end




%% ============================================================
% Output directory
% ============================================================

outputRoot = fullfile( ...
    analysisRoot, ...
    'SubjectResults', ...
    groupName, ...
    subj);

if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end


%% ========================================================================
% Save file 1
%
% Second-view inputs for later LG / GL FWCCA CV
% =========================================================================
Second_view_inputs_outputFile = fullfile( ...
    outputRoot, ...
    sprintf('%s_SecondViewInputs.mat', subj));

save( ...
     Second_view_inputs_outputFile, ...
     'AllTrainInputs', ...
    '-v7.3');

fprintf('\n');
fprintf('============================================================\n');
fprintf('Saving Second-view inputs for later LG / GL FWCCA CV \n');
fprintf('============================================================\n');
fprintf('%s\n', Second_view_inputs_outputFile);

%% ========================================================================
% Save file 2
%
% CCA / G-FWCCA training results
% =========================================================================

trainingOutputFile = fullfile( ...
    outputRoot, ...
    sprintf( ...
        '%s_CCA_GFWCCA_TrainingResults.mat', ...
        subj));


save( ...
    trainingOutputFile, ...
    'CCA_GFWCCA_TrainingResults', ...
    '-v7.3');

fprintf('\n');
fprintf('============================================================\n');
fprintf('Saving CCA / G-FWCCA training results\n');
fprintf('============================================================\n');
fprintf('%s\n', trainingOutputFile);

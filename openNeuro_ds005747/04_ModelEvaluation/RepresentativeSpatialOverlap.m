% =========================================================================
% RepresentativeSpatialOverlap.m
%
% Purpose:
%   Perform spatial analysis of the representative subjects identified
%   from the component-wise held-out correlation results.
%
%   For each selected ROI/component:
%     1. Extract the held-out seed-correlation and canonical maps.
%     2. Sign-align and normalize the maps to [0,1].
%     3. Compute voxel-wise distance to the identity line y = x.
%     4. Compute LG-FWCCA and GL-FWCCA improvements relative to CCA.
%     5. Retain the top 20% of positive improvements for each local method.
%     6. Identify common LG/GL improvement voxels and compute Jaccard overlap.
%     7. Determine the median spatial location of the common-improvement region.
%     8. Plot the held-out seed-correlation map and common-improvement region
%        in sagittal, coronal, and axial views.
%
% ROIs/components:
%   AnG  : component 4
%   PCgG : component 2
%   PCu  : component 1
%
% Inputs:
%   ../03_ModelTesting/Results/rsfMRI_TestingResults.mat
%   Results/RepresentativeSubjectResults.mat
%
% Outputs:
%   Results/RepresentativeSpatialOverlap/SpatialOverlapResults.mat
%   3 ROIs x 2 map types x 3 views = 18 figures
%
% =========================================================================

clear;
clc;
close all;


%% ========================================================================
% Project paths
% =========================================================================

% openNeuro_ds005747/04_ModelEvaluation
scriptDir = fileparts(mfilename('fullpath'));
datasetRoot = fileparts(scriptDir);
projectRoot = fileparts(datasetRoot);

spmDir = fullfile(projectRoot, 'spm_25');
resultsDir = fullfile(scriptDir, 'Results');

testFile = fullfile( ...
    datasetRoot, ...
    '03_ModelTesting', ...
    'Results', ...
    'rsfMRI_TestingResults.mat');

repFile = fullfile( ...
    resultsDir, ...
    'RepresentativeSubjectResults.mat');

outputDir = fullfile( ...
    resultsDir, ...
    'RepresentativeSpatialOverlap');

pdfDir = fullfile(outputDir, 'pdf');
pngDir = fullfile(outputDir, 'png');
figDir = fullfile(outputDir, 'fig');

if ~exist(outputDir,'dir'), mkdir(outputDir); end
if ~exist(pdfDir,'dir'), mkdir(pdfDir); end
if ~exist(pngDir,'dir'), mkdir(pngDir); end
if ~exist(figDir,'dir'), mkdir(figDir); end


%% ========================================================================
% SPM setup
% =========================================================================

addpath(spmDir);
spm('defaults','fmri');
spm_jobman('initcfg');


%% ========================================================================
% Load held-out maps and representative-subject results
% =========================================================================

S_test = load(testFile, 'AllTestMaps');
AllTestMaps = S_test.AllTestMaps;

S_rep = load(repFile, 'RepresentativeSubjectResults');
RepresentativeSubjectResults = S_rep.RepresentativeSubjectResults;


%% ========================================================================
% Neuromorphometrics atlas
% =========================================================================

atlas_file = fullfile( ...
    spmDir, 'atlas', 'Neuromorphometrics', ...
    'Neuromorphometrics.nii');

atlas_vol = spm_vol(atlas_file);
atlas_data = spm_read_vols(atlas_vol);

vol_size_3d = atlas_vol.dim;
vol_size = prod(vol_size_3d);


%% ========================================================================
% Load and reslice SPM T1 anatomical template to atlas grid
% =========================================================================

t1_file = fullfile(spmDir, 'canonical', 'avg152T1.nii');
t1_vol = spm_vol(t1_file);

same_dim = isequal(t1_vol.dim, atlas_vol.dim);
same_affine = max(abs(t1_vol.mat(:)-atlas_vol.mat(:))) < 1e-5;

if ~(same_dim && same_affine)

    flags = struct( ...
        'interp',4, ...
        'wrap',[0 0 0], ...
        'mask',0, ...
        'which',1, ...
        'mean',0, ...
        'prefix','atlas_');

    spm_reslice(char(atlas_file,t1_file), flags);

    [t1_dir,t1_name,t1_ext] = fileparts(t1_file);

    t1_resliced_file = fullfile( ...
        t1_dir, ['atlas_' t1_name t1_ext]);

else
    t1_resliced_file = t1_file;
end

t1_data = spm_read_vols(spm_vol(t1_resliced_file));


%% ========================================================================
% ROI settings
% =========================================================================

roi_list = {'AnG','PCgG','PCu'};

roi_labels = { ...
    [106 107], ...       % AnG
    [166 167], ...       % PCgG
    [168 169]};          % PCu

topPercent = 20;

SpatialOverlapResults = struct();


%% ========================================================================
% ROI loop
% =========================================================================

for roi_idx = 1:numel(roi_list)

    roiName = roi_list{roi_idx};

    subjectName = RepresentativeSubjectResults.(roiName).subject;
    k = RepresentativeSubjectResults.(roiName).component;

    fprintf('\n============================================================\n');
    fprintf('%s | %s | k = %d\n', roiName, subjectName, k);
    fprintf('============================================================\n');


    %% --------------------------------------------------------------------
    % Representative subject
    % ---------------------------------------------------------------------

    roiResults = AllTestMaps.(roiName);
    subjectNames = {roiResults.subject};
    sub_id = find(strcmp(subjectNames,subjectName),1);

    if isempty(sub_id)
        error('Cannot find %s for ROI %s.', subjectName, roiName);
    end


    %% --------------------------------------------------------------------
    % Extract held-out maps
    % ---------------------------------------------------------------------

    seed    = roiResults(sub_id).rseedtest_vals;
    CCA_map = roiResults(sub_id).CCA(:,k);
    LG_map  = roiResults(sub_id).LG_FWCCA(:,k);
    GL_map  = roiResults(sub_id).GL_FWCCA(:,k);


    %% --------------------------------------------------------------------
    % Sign-align canonical maps with held-out seed map
    % ---------------------------------------------------------------------

    if corr(seed,CCA_map,'Type','Pearson') < 0, CCA_map = -CCA_map; end
    if corr(seed,LG_map,'Type','Pearson')  < 0, LG_map  = -LG_map;  end
    if corr(seed,GL_map,'Type','Pearson')  < 0, GL_map  = -GL_map;  end


    %% --------------------------------------------------------------------
    % Normalize to [0,1]
    % ---------------------------------------------------------------------

    seed_norm = normalize_to_unit_interval(seed);
    CCA_norm  = normalize_to_unit_interval(CCA_map);
    LG_norm   = normalize_to_unit_interval(LG_map);
    GL_norm   = normalize_to_unit_interval(GL_map);


    %% --------------------------------------------------------------------
    % Distance to identity line and improvement relative to CCA
    % ---------------------------------------------------------------------

    dist_CCA = abs(CCA_norm-seed_norm)/sqrt(2);
    dist_LG  = abs(LG_norm-seed_norm)/sqrt(2);
    dist_GL  = abs(GL_norm-seed_norm)/sqrt(2);

    gain_LG = dist_CCA-dist_LG;
    gain_GL = dist_CCA-dist_GL;


    %% --------------------------------------------------------------------
    % Top 20% among positive improvements
    % ---------------------------------------------------------------------

    positive_LG = gain_LG(gain_LG > 0);
    positive_GL = gain_GL(gain_GL > 0);

    if isempty(positive_LG) || isempty(positive_GL)
        warning('%s has no positive LG or GL improvement.', roiName);
        continue;
    end

    thresholdPercentile = 100-topPercent;

    thr_LG = prctile(positive_LG, thresholdPercentile);
    thr_GL = prctile(positive_GL, thresholdPercentile);

    maskImprove_LG = gain_LG >= thr_LG;
    maskImprove_GL = gain_GL >= thr_GL;

    overlapMask = maskImprove_LG & maskImprove_GL;
    unionMask   = maskImprove_LG | maskImprove_GL;

    nLG      = nnz(maskImprove_LG);
    nGL      = nnz(maskImprove_GL);
    nOverlap = nnz(overlapMask);
    nUnion   = nnz(unionMask);

    if nUnion > 0
        jaccard = nOverlap/nUnion;
    else
        jaccard = NaN;
    end

    fprintf('LG top positive voxels = %d\n', nLG);
    fprintf('GL top positive voxels = %d\n', nGL);
    fprintf('Overlap voxels         = %d\n', nOverlap);
    fprintf('Jaccard overlap        = %.4f\n', jaccard);


    %% --------------------------------------------------------------------
    % Map ROI vectors into atlas space
    % ---------------------------------------------------------------------

    roiMask3D = ismember(atlas_data, roi_labels{roi_idx});
    roiIndices = find(roiMask3D);

    if numel(roiIndices) ~= numel(seed)
        error('ROI voxel mismatch for %s.', roiName);
    end

    seedVol = zeros(vol_size,1);
    seedVol(roiIndices) = seed_norm;
    seed3D = reshape(seedVol,vol_size_3d);

    overlapVol = zeros(vol_size,1);
    overlapVol(roiIndices) = double(overlapMask);
    overlap3D = reshape(overlapVol,vol_size_3d);


    %% --------------------------------------------------------------------
    % Median spatial location of common-improvement voxels
    % ---------------------------------------------------------------------

    overlap_indices = roiIndices(overlapMask);

    if isempty(overlap_indices)
        warning('%s has no LG/GL overlap voxels.', roiName);
        continue;
    end

    [xv,yv,zv] = ind2sub(vol_size_3d,overlap_indices);

    x_center = round(median(xv));
    y_center = round(median(yv));
    z_center = round(median(zv));

    world_center = atlas_vol.mat * ...
        [x_center; y_center; z_center; 1];

    world_center = world_center(1:3)';

    fprintf('Center voxel = (%d,%d,%d)\n', ...
        x_center,y_center,z_center);

    fprintf('Center MNI   = (%.1f, %.1f, %.1f)\n', ...
        world_center);


    %% --------------------------------------------------------------------
    % Store spatial-overlap results
    % ---------------------------------------------------------------------

    SpatialOverlapResults.(roiName).subject = subjectName;
    SpatialOverlapResults.(roiName).component = k;
    SpatialOverlapResults.(roiName).topPercent = topPercent;

    SpatialOverlapResults.(roiName).thresholdLG = thr_LG;
    SpatialOverlapResults.(roiName).thresholdGL = thr_GL;

    SpatialOverlapResults.(roiName).distanceCCA = dist_CCA;
    SpatialOverlapResults.(roiName).distanceLG = dist_LG;
    SpatialOverlapResults.(roiName).distanceGL = dist_GL;

    SpatialOverlapResults.(roiName).gainLG = gain_LG;
    SpatialOverlapResults.(roiName).gainGL = gain_GL;

    SpatialOverlapResults.(roiName).maskLG = maskImprove_LG;
    SpatialOverlapResults.(roiName).maskGL = maskImprove_GL;
    SpatialOverlapResults.(roiName).overlapMask = overlapMask;
    SpatialOverlapResults.(roiName).overlapVoxelIndices = overlap_indices;

    SpatialOverlapResults.(roiName).nLG = nLG;
    SpatialOverlapResults.(roiName).nGL = nGL;
    SpatialOverlapResults.(roiName).nOverlap = nOverlap;
    SpatialOverlapResults.(roiName).nUnion = nUnion;
    SpatialOverlapResults.(roiName).jaccard = jaccard;

    SpatialOverlapResults.(roiName).centerVoxel = ...
        [x_center y_center z_center];

    SpatialOverlapResults.(roiName).centerMNI = world_center;

    SpatialOverlapResults.(roiName).seedNormalized = seed_norm;
    SpatialOverlapResults.(roiName).CCANormalized = CCA_norm;
    SpatialOverlapResults.(roiName).LGNormalized = LG_norm;
    SpatialOverlapResults.(roiName).GLNormalized = GL_norm;


    %% ====================================================================
    % Prepare sagittal, coronal, and axial slices
    % =====================================================================

    viewNames = {'Sagittal','Coronal','Axial'};

    viewCoords = { ...
        sprintf('x%d',x_center), ...
        sprintf('y%d',y_center), ...
        sprintf('z%d',z_center)};

    anatSlices = { ...
        squeeze(t1_data(x_center,:,:))', ...
        squeeze(t1_data(:,y_center,:))', ...
        squeeze(t1_data(:,:,z_center))'};

    seedSlices = { ...
        squeeze(seed3D(x_center,:,:))', ...
        squeeze(seed3D(:,y_center,:))', ...
        squeeze(seed3D(:,:,z_center))'};

    roiSlices = { ...
        squeeze(roiMask3D(x_center,:,:))', ...
        squeeze(roiMask3D(:,y_center,:))', ...
        squeeze(roiMask3D(:,:,z_center))'};

    overlapSlices = { ...
        squeeze(overlap3D(x_center,:,:))', ...
        squeeze(overlap3D(:,y_center,:))', ...
        squeeze(overlap3D(:,:,z_center))'};


    %% ====================================================================
    % Map type 1: held-out seed-correlation map
    % =====================================================================

    for v = 1:3

        fig = figure( ...
            'Color','w', ...
            'Units','inches', ...
            'Position',[1 1 4.5 4]);

        ax = axes(fig);

        anat = prepare_T1_slice(anatSlices{v});
        image(ax,repmat(anat,[1 1 3]));

        axis(ax,'image');
        axis(ax,'off');
        hold(ax,'on');

        hSeed = imagesc(ax,seedSlices{v},[0 1]);
        colormap(ax,hot);

        set(hSeed, ...
            'AlphaData', ...
            0.80*double(roiSlices{v}));

        axis(ax,'image');
        axis(ax,'off');

        if v == 3
            cb = colorbar(ax);
            cb.Ticks = [0 0.5 1];
            cb.Label.String = 'Normalized correlation';
            cb.FontSize = 18;
            cb.Label.FontSize = 18;
        end

        baseName = sprintf( ...
            '%s_%s_SeedMap_%s_%s', ...
            roiName,subjectName,viewNames{v},viewCoords{v});

        exportgraphics(fig, ...
            fullfile(pdfDir,[baseName '.pdf']), ...
            'ContentType','vector');

        exportgraphics(fig, ...
            fullfile(pngDir,[baseName '.png']), ...
            'Resolution',300);

        savefig(fig,fullfile(figDir,[baseName '.fig']));

    end


    %% ====================================================================
    % Map type 2: common LG/GL positive-improvement region
    % =====================================================================

    for v = 1:3

        fig = figure( ...
            'Color','w', ...
            'Units','inches', ...
            'Position',[1 1 4 4]);

        ax = axes(fig);

        anat = prepare_T1_slice(anatSlices{v});
        image(ax,repmat(anat,[1 1 3]));

        axis(ax,'image');
        axis(ax,'off');
        hold(ax,'on');

        contour( ...
            ax, ...
            double(overlapSlices{v}), ...
            [0.5 0.5], ...
            'k-', ...
            'LineWidth',2.0);

        axis(ax,'image');
        axis(ax,'off');

        baseName = sprintf( ...
            '%s_%s_k%d_CommonImprovement_%s_%s', ...
            roiName,subjectName,k,viewNames{v},viewCoords{v});

        exportgraphics(fig, ...
            fullfile(pdfDir,[baseName '.pdf']), ...
            'ContentType','vector');

        exportgraphics(fig, ...
            fullfile(pngDir,[baseName '.png']), ...
            'Resolution',300);

        savefig(fig,fullfile(figDir,[baseName '.fig']));

    end

end


%% ========================================================================
% Save numerical spatial-overlap results
% =========================================================================

saveFile = fullfile( ...
    outputDir, ...
    'SpatialOverlapResults.mat');

save(saveFile, ...
    'SpatialOverlapResults', ...
    '-v7.3');

fprintf('\nSpatial-overlap results saved to:\n%s\n',saveFile);


%% ========================================================================
% Helper function: normalize T1 slice for display
% =========================================================================

function anat = prepare_T1_slice(anat)

    anat = double(anat);

    validVals = anat(isfinite(anat) & anat > 0);

    if isempty(validVals)
        anat = zeros(size(anat));
        return;
    end

    lowVal  = prctile(validVals,1);
    highVal = prctile(validVals,99);

    anat = (anat-lowVal) ./ max(highVal-lowVal,eps);
    anat = min(max(anat,0),1);

    % Slightly lighten anatomical background
    anat = 0.20 + 0.80*anat;

end

function [final_mask, idx_mask, idx_non, info] = generate_clean_mask( ...
    mask_file, rc1_file, rc2_file, rc3_file, opts)

%   GENERATE_CLEAN_MASK Construct a clean brain mask for CCA analysis.
%
%   [final_mask, idx_mask, idx_non, info] = GENERATE_CLEAN_MASK(
%       mask_file, rc1_file, rc2_file, rc3_file, opts)
%   constructs a voxel-level brain mask suitable for spatial CCA / FWCCA
%   analysis by retaining gray matter (and optionally white matter), while
%   excluding CSF, ventricles, the brain stem, and other non-neuronal regions.
%
%   INPUTS:
%     mask_file : SPM-generated brain mask (e.g., mask.nii).
%     rc1_file  : Resliced gray-matter probability map (c1*.nii).
%     rc2_file  : Resliced white-matter probability map (c2*.nii).
%     rc3_file  : Resliced CSF probability map (c3*.nii).
%     opts      : Optional structure specifying thresholds and exclusion rules.
%
%   OUTPUTS:
%     final_mask : 3D logical mask defining voxels retained for analysis.
%     idx_mask   : Linear indices of voxels inside final_mask.
%     idx_non    : Linear indices of voxels outside final_mask.
%     info       : Structure containing intermediate masks and configuration options.
%
%   NOTES:
%     - Tissue probability thresholds are applied to GM, WM, and CSF maps.
%     - Anatomical exclusion is performed using the Neuromorphometrics atlas.
%     - Optional removal of small connected components is supported.
%
%   This function defines the spatial domain (voxels) used for FWCCA analysis.

% ---------- Default parameters ----------
if nargin < 5, opts = struct(); end
opts = set_default(opts, 'gm_thr', 0.60);
opts = set_default(opts, 'wm_thr', 0.60);
opts = set_default(opts, 'use_wm', false);
opts = set_default(opts, 'csf_thr', 0.10);
opts = set_default(opts, 'exclude_labels', {'Brain Stem','4th Ventricle','3rd Ventricle', ...
     'Left Lateral Ventricle','Right Lateral Ventricle'});
opts = set_default(opts, 'atlas_name', 'Neuromorphometrics');
opts = set_default(opts, 'min_cluster', 0);


% ---------- Load input volumes ----------
Vmask = spm_vol(mask_file);
refMask = spm_read_vols(Vmask) > 0;      % initial brain mask
rc1 = spm_read_vols(spm_vol(rc1_file));  % GM
rc2 = spm_read_vols(spm_vol(rc2_file));  % WM
rc3 = spm_read_vols(spm_vol(rc3_file));  % CSF

% ---------- Tissue-based masking ----------
gm_mask = rc1 > opts.gm_thr;
wm_mask = opts.use_wm * (rc2 > opts.wm_thr);
csf_mask = rc3 > opts.csf_thr;
brain_mask = (gm_mask | wm_mask) & ~csf_mask;
base_mask = refMask & brain_mask;

% ---------- Atlas-based exclusion ----------
atlas_excl = false(size(base_mask));
try
    A = spm_atlas('load', opts.atlas_name);
    label_names = {A.labels.name};
    for i = 1:numel(opts.exclude_labels)
        label = opts.exclude_labels{i};
        label_idx = find(strcmp(label_names, label));
        if isempty(label_idx)
            warning('[!] Label not found: "%s"', label);
            continue;
        end

        
        label_val = A.labels(label_idx).index;
        Vm = spm_atlas('mask', A, label);
        if isstruct(Vm) && strcmp(Vm.fname, 'Neuromorphometrics_mask.nii')
            Vm.fname = fullfile(spm('Dir'), 'atlas', 'Neuromorphometrics', 'Neuromorphometrics.nii');
        end

        [p, n, e] = fileparts(Vm.fname);
        rfile = fullfile(p, ['r' n e]);

        if ~isfile(rfile)
            spm_reslice(char(Vmask.fname, Vm.fname), struct('which',1,'interp',0,'mean',false));
        end

        Vr = spm_vol(rfile);
        Y  = spm_read_vols(Vr);
        mLbl = (round(Y) == label_val);

        atlas_excl = atlas_excl | mLbl;
        fprintf('[✓] Excluded label "%s": %d voxels masked.\n', label, nnz(mLbl));
    end
catch ME
    warning('[!] Atlas exclusion skipped: %s', ME.message);
end

% ---------- Final mask ----------
mask_no_atlas = base_mask & ~atlas_excl;
final_mask = mask_no_atlas; %3d logical mask


% ---------- Remove small clusters (optional) ----------
if opts.min_cluster > 0
    try
        CC = bwconncomp(final_mask, 6);
        sizes = cellfun(@numel, CC.PixelIdxList);
        small = sizes < opts.min_cluster;
        final_mask(cat(1, CC.PixelIdxList{small})) = false;
    catch
        warning('[!] Skipped cluster removal (Image Toolbox required).');
    end
end

% ---------- Index sets ----------
idx_mask = find(final_mask);               % linear indice in MIN150
idx_non  = setdiff(1:numel(final_mask), idx_mask);

% ---------- Diagnostic information ----------info = struct();
info.opts          = opts;
info.refMask       = refMask;
info.gm_mask       = gm_mask;
info.wm_mask       = wm_mask;
info.csf_mask      = csf_mask;
info.brain_mask    = brain_mask;
info.base_mask     = base_mask;
info.atlas_excl    = atlas_excl;
info.mask_no_atlas = mask_no_atlas;

end

function s = set_default(s, field, val)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = val;
end
end
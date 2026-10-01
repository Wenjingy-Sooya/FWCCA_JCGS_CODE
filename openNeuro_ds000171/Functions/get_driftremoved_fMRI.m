function Y_filtered = get_driftremoved_fMRI(spm_mat_path, nii_path)

%GET_DRIFTREMOVED_FMRI Remove low-frequency temporal drift from fMRI time series.
%
%
%   INPUTS:
%     spm_mat_path : path to SPM.mat containing the design matrix and filter
%                    settings (SPM.xX.K).
%     nii_path     : path to a 4D preprocessed fMRI NIfTI file (e.g., swar*.nii).
%
%   OUTPUT:
%     Y_filtered   : drift-removed fMRI data in matrix form
%                    [nVoxels x nTimePoints].
%
%   NOTES:
%     - Temporal filtering is performed using SPM's spm_filter function,
%       ensuring consistency with standard SPM preprocessing.
%     - No spatial smoothing or whitening is applied in this function.
%
%   This function is used to construct the first view (X) in FWCCA analyses.
  
% Load design matrix and filter struct
    load(spm_mat_path, 'SPM');

    % --- Step 1: Read 4D fMRI image ---
    V = spm_vol(nii_path);                      % Get header for all volumes
    Y_raw = spm_read_vols(V);                   % 4D: [x y z t]
    Y_2D = reshape(Y_raw, [], size(Y_raw, 4));  % [nVoxels x nTime]

    % --- Step 2: High-pass filter (drift removal) ---
    Y_filtered = spm_filter(SPM.xX.K, Y_2D')';  % [nVox x nTime]


end

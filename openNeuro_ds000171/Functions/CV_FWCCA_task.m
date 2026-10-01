
function [best_hAndr, mean_count] = CV_FWCCA_task(seedvoxel_indices, ...
                                                  X1_seed, X2_seed, vol_3d, ...
                                                  atlas, nii_path, Wg, ...
                                                  Kcomp, h_list, radius_list, ...
                                                  num_fold, TR, mode, ridgeTol, svdTol)
 
% Cross-validate the bandwidth h and neighborhood radius r for LG-FWCCA and
% GL-FWCCA using K-fold cross-validation over voxels.
%
% INPUTS:
%   seedvoxel_indices : [N x 1] linear voxel indices defining the seed or
%                       pseudo-labeled region used for evaluation.
%   X1_seed, X2_seed  : [N x T] data matrices (voxels × time), restricted to
%                       the seed voxels for the two CCA views.
%   vol_3d            : [1 x 3] 3D volume size (e.g., [91 109 91]), used for
%                       index-to-coordinate conversion.
%   atlas             : Atlas structure loaded via spm_atlas
%                       (e.g., spm_atlas('load','Neuromorphometrics')).
%   nii_path          : Path to a representative NIfTI file, used to obtain
%                       voxel-to-MNI coordinate transforms.
%   Wg                : [T x T] initial global temporal weight matrix.
%   Kcomp             : Number of canonical components to estimate.
%   h_list            : List of candidate bandwidth values.
%   radius_list       : List of candidate neighborhood radius values.
%   num_fold          : Number of cross-validation folds.
%   TR                : Repetition time (in seconds).
%   mode              : FWCCA mixing mode:
%                       'prescale'  (global → local, GL-FWCCA) or
%                       'postscale' (local → global, LG-FWCCA).
%
% OUTPUTS:
%   mean_count : [num_h x num_r x Kcomp] average number of distinct anatomical
%                regions identified across folds, for each (h, r) pair.
%   best_hAndr : [1 x 2] selected (h, r) pair that maximizes the mean region
%                count aggregated across canonical components.

% --------------------------------------------------------------
% Initialization
% --------------------------------------------------------------
num_h = length(h_list);
num_r = length(radius_list);
[V1, ~] = size(X1_seed);
vol_size = vol_3d;                 % e.g., [91 109 91]
numofvol = prod(vol_size);         % total voxel count
V = spm_vol(nii_path);             % For MNI transform

fprintf('=== Starting CV for task-fMRI FWCCA ===\n');
fprintf('Total voxels: %d | h_list: %d | r_list: %d | folds: %d\n', ...
        V1, num_h, num_r, num_fold);

% --------------------------------------------------------------
% Step 1: Create CV folds
% --------------------------------------------------------------
fold_sizes = repmat(floor(V1 / num_fold), num_fold, 1);
fold_sizes(1:mod(V1, num_fold)) = fold_sizes(1:mod(V1, num_fold)) + 1;

fold_ids = zeros(V1, 1);
current_idx = 1;
for f = 1:num_fold
    fold_range = current_idx:(current_idx + fold_sizes(f) - 1);
    fold_ids(fold_range) = f;
    current_idx = current_idx + fold_sizes(f);
end

Exclude_idx = false(V1, num_fold);
for fold = 1:num_fold
    Exclude_idx(:, fold) = (fold_ids == fold);
end

% Allocate memory
NumofRegion_Cube = zeros(num_h, num_r, num_fold, Kcomp);

% --------------------------------------------------------------
% Step 2: Loop over h and r
% --------------------------------------------------------------
for h_idx = 1:num_h
    h = h_list(h_idx);
    for r_idx = 1:num_r
        r_current = radius_list(r_idx);
        fprintf('--- h = %.2f | r = %.2f ---\n', h, r_current);

        % Precompute spatial kernel (based on radius & h)
        K = makeLocalKernel(size(X1_seed, 2), r_current, 'gaussian', h, TR);

        for CV_idx = 1:num_fold
            fprintf('   Fold %d/%d ... ', CV_idx, num_fold);
            mask_out = Exclude_idx(:, CV_idx);

            % Split train/test
            X1_train = X1_seed(~mask_out, :);
            X2_train = X2_seed(~mask_out, :);
            X1_test  = X1_seed(mask_out, :);
            CV_test_indices = seedvoxel_indices(mask_out);

            % --- Train FWCCA on training data ---
            res = Run_FWCCA_Family(X1_train, X2_train, Wg, Wg, K, K, Kcomp, 100, 1e-5, mode, ridgeTol, svdTol);

            % --- Project held-out data ---
            maps_test = (X1_test - res.X1_localMean) * res.Ak_raw;  % [V_seed x Kcomp]

            for id_k = 1:Kcomp
                comp = abs(maps_test(:, id_k));
                comp_mask = comp / max(comp);

                % Sort and find top voxels
                [~, sort_idx] = sort(comp_mask, 'descend');
                all_region_names = {};
                num_top_pick = floor(0.0005 * numel(sort_idx));
                % --- Extract top regions (top 0.05% voxels) ---
                for k = 1:num_top_pick
                    max_voxel_idx = CV_test_indices(sort_idx(k));
                    [x, y, z] = ind2sub(vol_size, max_voxel_idx);
                    mni_coords = V(1).mat * [x; y; z; 1];
                    xY = struct('xyz', mni_coords(1:3), 'def', 'sphere', 'spec', 1);
                    region = spm_atlas('query', atlas, xY);
                    all_region_names{end+1} = strjoin(region, ', ');
                end

                % --- Clean names and count unique regions ---
                flat_names = strsplit(strjoin(all_region_names, ', '), ', ');
                flat_names = strtrim(flat_names);
                flat_names = flat_names(~cellfun(@isempty, flat_names));
                flat_names = flat_names(~contains(flat_names, 'Unknown'));

                NumofRegion_Cube(h_idx, r_idx, CV_idx, id_k) = numel(unique(flat_names));
            end
            fprintf('done.\n');
        end
    end
end

%% --------------------------------------------------------------
% Step 3: Average across folds
%
% Preserve dimensions explicitly:
%   [num_h x num_r x Kcomp]
% --------------------------------------------------------------

mean_count = mean(NumofRegion_Cube, 3);

mean_count = reshape( ...
    mean_count, ...
    num_h, ...
    num_r, ...
    Kcomp);

sum_K_count = sum(mean_count, 3);   % [num_h x num_r]

if isscalar(sum_K_count)
    % Only one (h,r) tested
    max_count = sum_K_count;
    best_hAndr = [h_list(1), radius_list(1)];
else
    % Multiple h and/or r
    [max_count, max_idx] = max(sum_K_count(:));
    [best_h_idx, best_r_idx] = ind2sub(size(sum_K_count), max_idx);

    % Protect against scalar h_list or r_list
    if isscalar(h_list)
        best_h_idx = 1; 
    end

    if isscalar(radius_list)
        best_r_idx = 1;
    end

    best_hAndr = [h_list(best_h_idx), radius_list(best_r_idx)];
end
fprintf('=== Best parameters found: h = %.3f | r = %.3f (mean count = %.2f) ===\n', ...
        best_hAndr(1), best_hAndr(2), max_count);
end

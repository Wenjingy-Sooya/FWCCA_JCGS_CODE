
function [X2, V_active_pattern] = BuildActiveCCAViews(Y_4D, mask_3D, Z0_3D)
%BUILDACTIVECCAVIEWS Construct a hybrid second view for spatial CCA.
%
%   [X2, V_active_pattern] = BuildActiveCCAViews(Y_4D, mask_3D, Z0_3D)
%   constructs a hybrid 4D representation in which voxels deemed active
%   (based on a pseudo-label map) are replaced by locally averaged time
%   series from their spatial neighbors, while inactive voxels are assigned
%   a global mean signal.
%
%   INPUTS:
%     Y_4D    : [nx × ny × nz × T] fMRI data array.
%     mask_3D : [nx × ny × nz] binary brain mask.
%     Z0_3D   : [nx × ny × nz] binary activation indicator (e.g., derived
%               from GLM-based pseudo-labels).
%
%   OUTPUTS:
%     X2              : [nx × ny × nz × T] hybrid-smoothed 4D volume defining
%                       the second CCA view.
%     V_active_pattern: [nx × ny × nz] logical mask indicating voxels treated
%                       as locally active.
%
%   NOTES:
%     - Active voxels are first dilated using a 3×3×3 structuring element to
%       define a local activation neighborhood.
%     - For active voxels, the time series is replaced by the average of
%       its 4-connected in-plane neighbors.
%     - For inactive voxels, a global mean time series is used.
%
%   This function implements the second-view construction strategy used in
%   the FWCCA framework for spatial fMRI analysis.


[nx, ny, nz, nTime] = size(Y_4D);

% --- Step 1: Get active voxel region
se = strel('cube', 3);
Z0_dilated = imdilate(Z0_3D, se);     % 3x3x3 neighborhood
V_active_pattern = (Z0_dilated) & mask_3D;

[x_act, y_act, z_act] = ind2sub([nx, ny, nz], find(V_active_pattern));

% --- Step 2: Initialize X2 as global mean signal for all voxels
X2 = zeros(nx, ny, nz, nTime);

% Compute global mean signal over all voxels in the brain mask
global_signal = squeeze(mean(Y_4D .* mask_3D, [1 2 3], 'omitnan'));

for t = 1:nTime
    X2(:,:,:,t) = global_signal(t);  % Broadcast to full volume
end

% --- Step 3: Overwrite active voxels with 4-neighbor average
for i = 1:length(x_act)
    x = x_act(i); y = y_act(i); z = z_act(i);

    % Check bounds to avoid edge errors
    if x <= 1 || x >= nx || y <= 1 || y >= ny
        continue;  % Skip boundary voxels
    end

    n_sum = zeros(1, nTime);
    n_valid = 0;

    % 4-connected neighbors: left, right, front, back
    offsets = [ -1  0  0;
                 1  0  0;
                 0 -1  0;
                 0  1  0 ];

    for j = 1:size(offsets, 1)
        xn = x + offsets(j,1);
        yn = y + offsets(j,2);
        zn = z + offsets(j,3);

        if xn >= 1 && xn <= nx && ...
           yn >= 1 && yn <= ny && ...
           zn >= 1 && zn <= nz && ...
           mask_3D(xn, yn, zn)

            n_sum = n_sum + squeeze(Y_4D(xn, yn, zn, :))';
            n_valid = n_valid + 1;
        end
    end

    if n_valid > 0
        X2(x, y, z, :) = n_sum / n_valid;
    else
        X2(x, y, z, :) = Y_4D(x, y, z, :);  % fallback
    end

end


end
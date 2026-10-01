
function [X2, V_active_pattern] = BuildActiveCCAViews(Y_4D, mask_3D, Z0_3D)

% BuildHybridX2: 
%   Constructs a 4D volume where activated voxels use spatial neighbors,
%   and non-activated voxels use the global average signal.
%
% INPUTS:
%   - Y_4D     : [nx × ny × nz × T] fMRI data
%   - mask_3D  : [nx × ny × nz] brain mask
%   - Z0_3D    : [nx × ny × nz] binary activation (e.g., from GLM)
%   - roi_3D   : [nx x ny x nz] intrested region (e.g., prior information)
% OUTPUT:
%   - X2       : same size as Y_4D, hybrid-smoothed 4D volume



[nx, ny, nz, nTime] = size(Y_4D);

% --- Step 1: Get active voxel region
se = strel('cube',3);
Z0_dilated = imdilate(Z0_3D, se);     % 3x3x3 neighborhood


%roi_dilated = imdilate(roi_3D, se);  % 3x3x3 neighborhood
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

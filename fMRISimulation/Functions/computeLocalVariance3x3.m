function localVar = computeLocalVariance3x3(spatial_map)
% spatial_map: [H x W] canonical spatial map (reshaped from vector form)
% localVar: scalar value representing average local variance

    % Pad the image to handle edges
    padded_map = padarray(spatial_map, [1 1], 'replicate');

    % Get size
    [H, W] = size(spatial_map);
    local_variances = zeros(H, W);

    % Compute local variance for each pixel
    for i = 2:H+1
        for j = 2:W+1
            patch = padded_map(i-1:i+1, j-1:j+1);  % 3x3 neighborhood
            local_variances(i-1,j-1) = var(patch(:), 1);  % variance of 9 values
        end
    end

    % Average local variance across all voxels
    localVar = mean(local_variances(:));

end
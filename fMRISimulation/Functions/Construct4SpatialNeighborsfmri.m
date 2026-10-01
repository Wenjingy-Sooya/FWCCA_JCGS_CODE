
function Y = Construct4SpatialNeighborsfmri(X, grid_size)
    % Constructs Y(t, i, j) as the sum of the four spatial neighbors of X(t, i, j)
    % 
    % Input:
    %   X: fMRI data (T × V), where each row is a timepoint, each column is a voxel
    %   grid_size: [M, N] specifying the spatial size of the brain slice (2D)
    % Output:
    %   Y: Spatially transformed data (T × V)

    [T, V] = size(X);  % Number of timepoints (T) and voxels (V)
    M = grid_size(1);  % Grid rows
    N = grid_size(2);  % Grid columns
  

    for t = 1:T  % Loop over all timepoints
        % Reshape 1D vector X(t, :) into 2D spatial grid (M × N)
        X_t_2D = reshape(X(t, :), M, N);
        
        % Sum valid neighbors (no boundary corrections)
        BottomMat = [X_t_2D(2:M, :); zeros(1, N)]; % Correct
        TopMat = [zeros(1, N); X_t_2D(1:M-1, :)]; % Correct
        RightMat=[X_t_2D(:, 2:N) zeros(M,1) ];  % Right neighbor
        LeftMat=[zeros(M,1)  X_t_2D(:, 1:N-1)]; % Left neighbor
         % Compute Y_t_2D using spatial neighbor sum
       Y_t_2D = BottomMat + TopMat + RightMat + LeftMat;

        % Reshape Y_t_2D back to 1D and store it in Y
        Y(t, :) = reshape(Y_t_2D, 1, V);
    end
end
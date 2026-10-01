function [rows, cols] = block_from_center(center, sz, grid_size)
%% ============================================================
% BLOCK_FROM_CENTER Construct a rectangular spatial block from its center
%
% This function determines the row and column indices of a rectangular
% spatial block given its nominal center and block size. If the block
% would extend beyond the spatial grid, its starting position is shifted
% so that the complete block remains within the grid boundaries.
%
% This function is used to construct the simulated spatial signal regions
% after applying random spatial jitter to their center locations.
%
% Input
% -----
% center : 1-by-2 vector
%   Nominal center of the spatial block:
%
%       center = [row, column]
%
% sz : 1-by-2 vector
%   Size of the rectangular block:
%
%       sz = [height, width]
%
% grid_size : 1-by-2 vector
%   Size of the complete spatial grid:
%
%       grid_size = [nRows, nCols]
%
% Output
% ------
% rows : vector
%   Row indices occupied by the spatial block.
%
% cols : vector
%   Column indices occupied by the spatial block.
%
% Example
% -------
%   center = [9, 9];
%   sz = [6, 6];
%   grid_size = [20, 20];
%
%   [rows, cols] = block_from_center(center, sz, grid_size);
%
%% ============================================================

r = center(1);
c = center(2);

h = sz(1);
w = sz(2);

% Half-widths used to determine the nominal upper-left corner
rh = floor((h - 1) / 2);
ch = floor((w - 1) / 2);

% Nominal starting row and column, constrained by the lower boundaries
r0 = max(1, r - rh);
c0 = max(1, c - ch);

% Shift the block if necessary to keep it inside the upper boundaries
r0 = min(r0, grid_size(1) - h + 1);
c0 = min(c0, grid_size(2) - w + 1);

% Row and column indices of the complete block
rows = r0:(r0 + h - 1);
cols = c0:(c0 + w - 1);

end
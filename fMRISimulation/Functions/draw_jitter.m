function d = draw_jitter(jitRange)
%% ============================================================
% DRAW_JITTER Generate a random two-dimensional spatial jitter
%
% This function randomly generates row and column displacements within
% user-specified integer ranges. It is used to perturb the locations of
% the simulated spatial signal regions across Monte Carlo replicates.
%
% Input
% -----
% jitRange : 2-by-2 matrix
%   Integer ranges for the row and column displacements:
%
%       jitRange = [row_min row_max;
%                   col_min col_max]
%
% Output
% ------
% d : 1-by-2 vector
%   Random spatial displacement:
%
%       d = [drow, dcol]
%
% Example
% -------
%   jitRange = [-2 2; -2 2];
%   d = draw_jitter(jitRange);
%
%% ============================================================

drow = randi(jitRange(1,2) - jitRange(1,1) + 1) ...
       + jitRange(1,1) - 1;

dcol = randi(jitRange(2,2) - jitRange(2,1) + 1) ...
       + jitRange(2,1) - 1;

d = [drow, dcol];

end
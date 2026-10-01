
function K = makeLocalKernel(T, radius, mode, sigma, TR)
% K: T x T temporal smoothing kernel, row-stochastic
% - Each row t contains weights to smooth timepoint t using nearby points.
% - Supports 'uniform' or 'gaussian' weights.
%
% Inputs:
%   T      : number of timepoints
%   radius : number of neighbors on each side (in index units)
%   mode   : 'uniform' or 'gaussian'
%   sigma  : std of Gaussian in seconds (only used if mode='gaussian')
%   TR     : sampling interval (repetition time) in seconds

    if nargin < 3, mode = 'uniform'; end
    if nargin < 4, sigma = 1.0; end
    if nargin < 5, TR = 1.0; end  % default 1s

    K = spalloc(T, T, (2*radius+1)*T);  % preallocate sparse T×T

    for t = 1:T
        idx = max(1, t-radius) : min(T, t+radius);  % neighbors in index
        switch mode
            case 'uniform'
                w = ones(numel(idx),1);
            case 'gaussian'
                d = abs(idx - t) * TR;  % convert to seconds!
                w = exp(-(d.^2) / (2*sigma^2));
            otherwise
                error('Unknown mode "%s". Must be "uniform" or "gaussian".', mode);
        end
        w = w / sum(w);  % row-normalize
        K(t, idx) = w;
    end
end
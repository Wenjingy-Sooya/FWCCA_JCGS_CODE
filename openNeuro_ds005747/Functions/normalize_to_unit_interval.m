function norm_vals = normalize_to_unit_interval(vals, mask)
%NORMALIZE_TO_UNIT_INTERVAL Normalize array to [0, 1] range.
%
%   norm_vals = normalize_to_unit_interval(vals)
%   norm_vals = normalize_to_unit_interval(vals, mask)
%
%   Inputs:
%     - vals: numeric array (vector, matrix, or volume)
%     - mask (optional): logical mask of same size as vals
%                        Only normalize values within the mask
%
%   Output:
%     - norm_vals: array of same size as vals, scaled to [0, 1] (only in mask)

    if nargin < 2
        mask = isfinite(vals);  % default: all finite entries
    end

    norm_vals = zeros(size(vals));  % initialize output

    % Extract values to normalize
    valid_vals = vals(mask);

    % Handle empty or constant values
    if isempty(valid_vals) || range(valid_vals) < 1e-6
        norm_vals(mask) = 0;  % or set to 0.5 if you prefer
        return;
    end

    % Normalize
    min_val = min(valid_vals);
    max_val = max(valid_vals);
    norm_vals(mask) = (valid_vals - min_val) / (max_val - min_val);
end
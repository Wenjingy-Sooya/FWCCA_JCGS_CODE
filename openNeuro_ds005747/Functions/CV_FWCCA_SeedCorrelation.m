function [meanCumCorr, BestParamTable, CV_seedCorr] = ...
        CV_FWCCA_SeedCorrelation( ...
        seedvoxel_indices, ...
        X1, X2, Wg, ...
        Kmaxcomp, ...
        h_list, radius_list, ...
        num_fold, TR, mode, ridgeTol, svdTol)

%CV_LEAVEONELABELOUT_FWCCA_VERSION3
% Select temporal bandwidth h and neighborhood radius r for
% LG-FWCCA / GL-FWCCA using num-fold cross-validation over voxels.
%
% The cross-validation criterion measures the correspondence between
% the held-out seed-correlation map and the held-out spatial canonical
% maps. For each candidate (h,r), the absolute Pearson correlations
% are accumulated over the first k canonical components:
%
%   C_k(h,r)= sum_{j=1}^k |corr(r_seed,test, s_j,test)|.
%
% The optimal (h,r) is selected separately for each retained number
% of components k = 1,...,Kmaxcomp.
%
%
% INPUTS
% -------------------------------------------------------------------------
% seedvoxel_indices : [N x 1]
%     Linear voxel indices defining the atlas-based ROI.
%
% X1, X2 : [N x T]
%     Training data matrices for the two CCA views.
%     Rows correspond to voxels and columns correspond to time points.
%
% Wg : [T x T]
%     Global temporal weighting matrix.
%
% Kmaxcomp : scalar
%     Maximum number of canonical components to estimate.
%
% h_list : vector
%     Candidate temporal bandwidth values.
%
% radius_list : vector
%     Candidate neighborhood-radius values.
%
% num_fold : scalar
%     Number of voxel-wise cross-validation folds.
%
% TR : scalar
%     Repetition time in seconds.
%
% mode : char/string
%     Mixing mode used by Run_FWCCA_Family:
%         'prescale'
%         'postscale'
%
%
% OUTPUTS
% -------------------------------------------------------------------------
% meanCumCorr : [num_h x num_r x Kmaxcomp]
%     Mean cumulative absolute seed-map/canonical-map correlation
%     across cross-validation folds.
%
%     meanCumCorr(h,r,k) corresponds to:
%
%       mean_fold sum_{j=1}^k ...
%           |corr(r_seed,test, s_j,test)|.
%
% BestParamTable : [Kmaxcomp x 4] table
%     Optimal tuning parameters for each retained number of canonical
%     components. The table contains:
%
%       K       - Number of retained canonical components.
%       Best_h  - Selected temporal bandwidth.
%       Best_r  - Selected neighborhood radius.
%       MeanCumulativeAbsCorrelation
%               - Maximum mean cross-validated cumulative absolute
%                 correlation between the held-out seed-correlation map
%                 and the first K held-out canonical maps.
% CV_seedCorr : cell array
%     Held-out seed-correlation vectors for diagnostic purposes.
%
% -------------------------------------------------------------------------


%% ============================================================
% Basic dimensions
%% ============================================================

num_h = numel(h_list);
num_r = numel(radius_list);

[Nvox, T] = size(X1);

if size(X2,1) ~= Nvox || size(X2,2) ~= T
    error('X1 and X2 must have identical dimensions.');
end

if numel(seedvoxel_indices) ~= Nvox
    error(['Length of seedvoxel_indices must match the number ' ...
           'of rows in X1/X2.']);
end

if size(Wg,1) ~= T || size(Wg,2) ~= T
    error('Wg must be a T x T temporal weight matrix.');
end


%% ============================================================
% Construct voxel-wise K-fold partitions
%
% Voxels are assigned sequentially to approximately equal folds.
%% ============================================================

fold_sizes = repmat(floor(Nvox / num_fold), num_fold, 1);
fold_sizes(1:mod(Nvox,num_fold)) =  fold_sizes(1:mod(Nvox,num_fold)) + 1;
fold_ids = zeros(Nvox,1);

current_idx = 1;

for f = 1:num_fold

    fold_range =  current_idx : (current_idx + fold_sizes(f) - 1);

    fold_ids(fold_range) = f;

    current_idx =  current_idx + fold_sizes(f);

end


%% ------------------------------------------------------------
% Exclusion masks
%% ------------------------------------------------------------

Exclude_idx = false(Nvox,num_fold);

for f = 1:num_fold
    Exclude_idx(:,f) =  (fold_ids == f);
end


%% ============================================================
% Result containers
%
% Dimensions:
%   h x r x fold x retained-component-number
%% ============================================================

corrCube = nan(num_h, num_r,  num_fold, Kmaxcomp);
CV_seedCorr = cell(num_fold,1);

%% ============================================================
% Cross-validation over candidate h and r
%% ============================================================

for h_idx = 1:num_h

    h = h_list(h_idx);

    for r_idx = 1:num_r

        radius_current = ...
            radius_list(r_idx);


        %% ----------------------------------------------------
        % Construct temporal local kernel
        %% ----------------------------------------------------
        K = makeLocalKernel( T, radius_current,  'gaussian', h,  TR);
        %% ====================================================
        % Cross-validation folds
        %% ====================================================

        for CV_idx = 1:num_fold

            mask_out =  Exclude_idx(:,CV_idx);


            %% ------------------------------------------------
            % Split voxels into training and held-out sets
            %% ------------------------------------------------
            X1_train = X1(~mask_out,:);
            X2_train = X2(~mask_out,:);
            X1_test = X1(mask_out,:);
            %% =================================================
            % Held-out seed-correlation map
            %
            % For each held-out voxel, compute its correlation
            % with the mean held-out ROI time series.
            %% =================================================

            seedtest_signal = mean(X1_test,1);
            seedtest_signal_z = (seedtest_signal - mean(seedtest_signal)) / std(seedtest_signal);
            X1_test_z = (X1_test - mean(X1_test,2)) ./ std(X1_test,0,2);
            r_test = X1_test_z * seedtest_signal_z' / (T - 1);
            CV_seedCorr{CV_idx} =  r_test;


            %% =================================================
            % Fit FWCCA on training voxels
            %% =================================================

            res = Run_FWCCA_Family( X1_train, X2_train, Wg, Wg, ...
                                    K,  K, Kmaxcomp, 100, 1e-5, mode, ridgeTol, svdTol);


            %% =================================================
            % Project held-out voxels
            %% =================================================

            maps_test =  (X1_test - res.X1_localMean)  * res.Ak_raw;

            r_test_norm = normalize_to_unit_interval(r_test);
         
            corrVec = arrayfun(@(k) corr( ...
                     r_test_norm, ...
                     normalize_to_unit_interval(maps_test(:,k)), ...
                    'Type', 'Pearson'), 1:Kmaxcomp);

            %% =================================================
            % Cumulative absolute correlation
            %
            % Element k represents the sum across components
            % 1,...,k.
            %% =================================================

            absCorr = abs(corrVec);
            cumulativeCorr = cumsum(absCorr);

            corrCube(h_idx, r_idx, CV_idx, :) = ...
                reshape(cumulativeCorr, ...
                    1,1,1,Kmaxcomp);  
 
            fprintf('h = %.2f | r = %.2f | fold = %d/%d\n',  h, radius_current, CV_idx,  num_fold);
            
        end

    end

end


%% ============================================================
% Average criterion across cross-validation folds
%
% Preserve dimensions explicitly:
%   [num_h x num_r x Kmaxcomp]
%% ============================================================

meanCumCorr = ...
    mean( ...
        corrCube, ...
        3, ...
        'omitnan');

meanCumCorr = ...
    reshape( ...
        meanCumCorr, ...
        num_h, ...
        num_r, ...
        Kmaxcomp);


%% ============================================================
% Select optimal (h,r) separately for each retained K
%% ============================================================

BestParamTable = table( ...
    (1:Kmaxcomp)', ...
    nan(Kmaxcomp,1), ...
    nan(Kmaxcomp,1), ...
    nan(Kmaxcomp,1), ...
    'VariableNames', { ...
        'K', ...
        'Best_h', ...
        'Best_r', ...
        'MeanCumulativeAbsCorrelation'});

for k = 1:Kmaxcomp

    criterion_k = ...
        meanCumCorr(:,:,k);

    %% --------------------------------------------------------
    % Ignore invalid entries
    %% --------------------------------------------------------

    criterion_for_max = ...
        criterion_k;

    criterion_for_max( ...
        ~isfinite(criterion_for_max)) = ...
        -Inf;

    [best_value, linear_idx] = ...
        max(criterion_for_max(:));

    if isinf(best_value)

        warning( ...
            'No valid (h,r) combination found for K = %d.', ...
            k);

        continue;

    end

    [best_h_idx,best_r_idx] = ...
        ind2sub( ...
            size(criterion_for_max), ...
            linear_idx);

    %% --------------------------------------------------------
    % Store selected parameters
    %% --------------------------------------------------------

    BestParamTable.Best_h(k) = h_list(best_h_idx);

    BestParamTable.Best_r(k) = radius_list(best_r_idx);

    BestParamTable.MeanCumulativeAbsCorrelation(k) =  best_value;

end

end
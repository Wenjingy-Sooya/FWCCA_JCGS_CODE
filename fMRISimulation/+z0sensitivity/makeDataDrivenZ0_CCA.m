function Z0 = makeDataDrivenZ0_CCA(X1c, X2c, activeFrac)

% Data-driven z0 from a preliminary unweighted CCA.
% X1c, X2c are N x T centered data matrices.
% activeFrac controls the proportion of selected voxels.

[N, ~] = size(X1c);

% Built-in MATLAB CCA
% A and B are canonical coefficient matrices.
% U and V are canonical variates / spatial maps.
[A, B, r, U, V] = canoncorr(X1c, X2c);

% Use the first canonical variate from the first view as preliminary map
score0 = abs(U(:,1));     % N x 1

% Select the top activeFrac proportion of voxels
nActive = max(1, round(activeFrac * N));

[~, idxSort] = sort(score0, 'descend');

Z0 = false(N,1);
Z0(idxSort(1:nActive)) = true;

Z0 = double(Z0);

end
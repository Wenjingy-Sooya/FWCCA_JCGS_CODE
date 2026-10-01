function Z0 = makeRandomMatchedZ0(Z_true)

Z_true = logical(Z_true(:));
N = numel(Z_true);
nActive = sum(Z_true);

Z0 = false(N,1);
idx = randperm(N, nActive);
Z0(idx) = true;

Z0 = double(Z0);

end
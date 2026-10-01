function w0 = makeInitialWeightsFromZ0(X1c, Z0, T)

Z0 = double(Z0(:));

% Center z0 for Pearson/cosine correlation
Z0c = Z0 - mean(Z0);

Xc = X1c - mean(X1c,1);

normX = sqrt(sum(Xc.^2,1))';   % T x 1
normZ = norm(Z0c);

corrVec = abs(Xc' * Z0c) ./ max(normX * normZ, eps);

if sum(corrVec) > eps
    w0 = T * corrVec / sum(corrVec);
else
    w0 = ones(T,1);
end

end
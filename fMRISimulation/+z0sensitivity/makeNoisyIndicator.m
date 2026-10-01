function Z0 = makeNoisyIndicator(Z_true, fpRate, fnRate)

Z_true = logical(Z_true(:));
N = numel(Z_true);

Z0 = Z_true;

activeIdx = find(Z_true);
inactiveIdx = find(~Z_true);

nFN = round(fnRate * numel(activeIdx));
nFP = round(fpRate * numel(inactiveIdx));

if nFN > 0
    fnIdx = activeIdx(randperm(numel(activeIdx), nFN));
    Z0(fnIdx) = false;
end

if nFP > 0
    fpIdx = inactiveIdx(randperm(numel(inactiveIdx), nFP));
    Z0(fpIdx) = true;
end

Z0 = double(Z0);

end
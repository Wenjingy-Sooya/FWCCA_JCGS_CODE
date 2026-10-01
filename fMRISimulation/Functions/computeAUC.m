function auc = computeAUC(scores, labels)

    scores = abs(scores(:));
    labels = labels(:);

    [~, sort_idx] = sort(scores,'descend');
    sorted_labels = labels(sort_idx);

    P = sum(sorted_labels == 1);
    N = sum(sorted_labels == 0);

    if P == 0 || N == 0
        warning('AUC undefined: only one class present. Returning NaN.');
        auc = NaN;
        return;
    end

    tps = cumsum(sorted_labels == 1);
    fps = cumsum(sorted_labels == 0);

    TPR = tps / P;
    FPR = fps / N;

    auc = trapz(FPR,TPR);
end

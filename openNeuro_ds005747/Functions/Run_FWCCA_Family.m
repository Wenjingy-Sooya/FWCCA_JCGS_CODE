% =========================================================================
% Run_FWCCA_Family.m
%
% Core implementation of the feature-weighted canonical correlation
% analysis (FWCCA) family used throughout the simulation, task-based fMRI,
% and resting-state fMRI analyses.
%
% The function estimates canonical directions under global feature
% weighting, local kernel weighting, or their combination. Global feature
% weights are updated iteratively using the association between the current
% canonical representation and the original features until convergence.
%
% DATA ORGANIZATION
%   X1, X2 : n x T data matrices
%            rows    = samples (e.g., voxels in the fMRI analyses)
%            columns = features (time points in the fMRI analyses)
%
% INPUTS
%   X1, X2      : Two input data views [n x T].
%   W1g_ini     : Initial global feature-weight matrix for X1 [T x T].
%   W2g_ini     : Initial global feature-weight matrix for X2 [T x T].
%   K1, K2      : Local feature-kernel matrices [T x T].
%   numOfCCA    : Number of canonical pairs to estimate.
%   maxIter     : Maximum number of alternating CCA updates.
%   tol         : Convergence tolerance.
%   mixOrder    : Order of global and local feature weighting:
%                   'prescale'  -> W^(1/2) * K  (GL-FWCCA)
%                   'postscale' -> K * W^(1/2)  (LG-FWCCA)
%   ridgeTol    : Ridge regularization parameter for covariance matrices.
%   svdTol      : Singular-value truncation tolerance.
%
% METHOD SPECIAL CASES
%   CCA:
%       W1g_ini = W2g_ini = I,  K1 = K2 = I
%
%   G-FWCCA:
%       global feature weights with K1 = K2 = I
%
%   LG-FWCCA:
%       local weighting followed by global feature weighting
%       mixOrder = 'postscale'
%
%   GL-FWCCA:
%       global feature weighting followed by local weighting
%       mixOrder = 'prescale'
%
% OUTPUT
%   results : Structure containing the final FWCCA fit, including
%             canonical directions, raw-time-axis directions, canonical
%             correlations, canonical component scores, updated global
%             feature weights, effective feature weights, and convergence
%             diagnostics.
%
% Important output fields include:
%   Ak, Bk                 : Canonical directions in the weighted space.
%   Ak_raw, Bk_raw         : Canonical directions mapped to the raw feature
%                            axis.
%   rhoVec                 : Canonical correlations.
%   comp1, comp2           : Canonical component scores for the two views.
%   Wg                     : Updated global feature-weight matrix.
%   effectiveWeight        : Effective feature weights induced by the
%                            combined global/local transformation.
%   X1_localMean,
%   X2_localMean           : Local means used for centering.
%   mixOrder               : Applied global/local weighting order.
%   outerIter              : Number of self-reinforcing outer iterations.
%   converged              : Indicator of convergence.
%   finalError             : Final change in effective feature weights.
%
% CONVERGENCE
%   The outer self-reinforcing iteration stops when the Euclidean change in
%   effective feature weights is no greater than tol, or when the maximum
%   number of outer iterations is reached.
%
% NOTES
%   - Local kernels are row-normalized internally.
%   - Covariance matrices are ridge-regularized and stabilized using
%     truncated SVD.
%   - Higher-order canonical pairs are obtained by sequential deflation.
%   - Additional convergence histories are retained in the output structure
%     for diagnostic purposes.
%
% This function is shared by all three reproducibility modules:
%   1. fMRISimulation
%   2. openNeuro_ds000171
%   3. openNeuro_ds005747
% =========================================================================

function results = Run_FWCCA_Family(X1, X2, W1g_ini, W2g_ini, K1, K2, numOfCCA, maxIter, tol, mixOrder, ridgeTol, svdTol)

    [~, T] = size(X1);

    %Initialize
    error = inf;
    outerIter = 0;
    maxOuterIter = 200;              % limit total iterations
    W1g_up = W1g_ini;
    W2g_up = W2g_ini;
    eff_pre = ones(T,1);             % initial effective weights (uniform)

    errorHistory = nan(maxOuterIter,1);
    Akraw_pre = [];
    rho_pre = [];
    AkrawHistory = nan(maxOuterIter,1);
    AkrawAngleHistory = nan(maxOuterIter, numOfCCA);
    rhoHistory = nan(maxOuterIter,numOfCCA);
    subspaceMaxAngleHistory = nan(maxOuterIter,1);
    subspaceMeanAngleHistory = nan(maxOuterIter,1);
    rhoGapHistory = nan(maxOuterIter,1);
   
    while error > tol && outerIter < maxOuterIter

        outerIter = outerIter + 1;

        % Save current effective weight for convergence check
        eff_old = eff_pre;

        % Run one self-reinforcing FWCCA update
        results_up = run_FWCCA_SelfReinforcing_globalLocalMix_famlily( ...
            X1, X2, W1g_up, W2g_up, K1, K2, numOfCCA, maxIter, tol , mixOrder, ridgeTol, svdTol);

        % Update weights and effective weight
        W1g_up = results_up.Wg;
        W2g_up = results_up.Wg;
        eff_pre = results_up.effectiveWeight;
        
        Akraw_cur = results_up.Ak_raw;
        rho_cur = results_up.rhoVec;

        if ~isempty(rho_pre)
          rhoHistory(outerIter,:) = rho_cur;
        end

        if numel(rho_cur) >=2

          rhoGapHistory(outerIter) = abs(rho_cur(1) - rho_cur(2));

        end
         
        rho_pre = rho_cur;

        if ~isempty(Akraw_pre)

          % Orthonormal bases for the two column spaces

          Qpre = orth(Akraw_pre);
          Qcur = orth(Akraw_cur);

        % Principal angles between the two canonical subspaces
          singularValues = svd(Qpre' * Qcur);
          singularValues = min(max(singularValues,0),1);
          principalAngles = acosd(singularValues);
          subspaceMaxAngleHistory(outerIter) = max(principalAngles);
          subspaceMeanAngleHistory(outerIter) = mean(principalAngles);

       
         % Align the sign of each current canonical direction
         % with the corresponding direction from the previous iteration
         Akraw_cur_aligned = Akraw_cur;

          for k = 1:numOfCCA

             if Akraw_pre(:,k)' * Akraw_cur(:,k) < 0

               Akraw_cur_aligned(:,k) = -Akraw_cur(:,k);

             end          

              % Canonical angle in degrees
               cosineValue = abs(Akraw_pre(:,k)' * Akraw_cur(:,k)) / (norm(Akraw_pre(:,k)) * norm(Akraw_cur(:,k)));
               % Protect against numerical values slightly outside [-1,1]
               cosineValue = min(max(cosineValue, -1), 1);
               AkrawAngleHistory(outerIter,k) = acosd(cosineValue);
          end

       % Frobenius norm of the change
          AkrawHistory(outerIter) =norm(Akraw_cur_aligned - Akraw_pre, 'fro');
       % Store the sign-aligned version for the next comparison
          Akraw_pre = Akraw_cur_aligned;

       else

          Akraw_pre = Akraw_cur;

       end

      % Compute change
        error = norm(eff_pre - eff_old);
        errorHistory(outerIter) = error;        
        %fprintf('Iter %d: ||w_new - w_old|| = %.5f\n', outerIter, error);
     end

    results = results_up;  % return the final result
                           % Convergence information for the self-reinforcing FWCCA

    results.outerIter = outerIter;
    results.converged = (error <= tol);
    results.finalError = error;
    results.svdTol = svdTol;
    results.errorHistory = errorHistory(1:outerIter);
    results.AkrawHistory = AkrawHistory(2:outerIter);
    results.AkrawAngleHistory = AkrawAngleHistory(2:outerIter,:);
    results.rhoHistory=rhoHistory(2:outerIter,:);
    results.subspaceMaxAngleHistory = subspaceMaxAngleHistory(2:outerIter);
    results.subspaceMeanAngleHistory = subspaceMeanAngleHistory(2:outerIter);
    results.rhoGapHistory = rhoGapHistory(1:outerIter);

end

%%%%%%%%%%%%Helpers###################

function results = run_FWCCA_SelfReinforcing_globalLocalMix_famlily(X1, X2, W1g_pre, W2g_pre, K1, K2, numOfCCA, maxIter, tol, mixOrder, ridgeTol, svdTol)

% Global + Local Weighted CCA 
% X1, X2: [n x T]
% W1g, W2g: [T x T] diagonal (global weights)
% K1, K2 : [T x T] row-stochastic local kernels
% mixOrder (optional): 'preScale' (W^1/2*K, default) or 'postScale' (K*W^1/2)

    if nargin < 10 || isempty(mixOrder)
        mixOrder = 'preScale';  % default: W^1/2 * K
    end

    assert(~isempty(svdTol), 'svdTol must be provided.');

    % ---------------- Safety & normalization ----------------
    [n, T] = size(X1);
    assert(all(size(X2) == [n, T]), 'X1 and X2 must be n x T.');
    W1g = diag(diag(W1g_pre)); 
    W2g = diag(diag(W2g_pre));
    K1 = normalize_rows(K1); 
    K2 = normalize_rows(K2);

    % ---------------- Local mean centering ------------------
    m1 = mean(X1, 1);  m2 = mean(X2, 1);
    m1_loc = m1 * K1;  m2_loc = m2 * K2;
    X1locentered = X1 - repmat(m1_loc, n, 1);
    X2locentered = X2 - repmat(m2_loc, n, 1);

    % ---------------- Build L1/L2 with chosen order ----------
    D1sqrt = diag(sqrt(diag(W1g)));
    D2sqrt = diag(sqrt(diag(W2g)));

    switch lower(mixOrder)
        case 'prescale'   % W^1/2 * K  (scale then mix) 
            L1 = D1sqrt * K1;
            L2 = D2sqrt * K2;
        case 'postscale'  % K * W^1/2  (mix then scale)
            L1 = K1 * D1sqrt;
            L2 = K2 * D2sqrt;
        otherwise
            error('mixOrder must be ''preScale'' or ''postScale''.');
    end

    % ---------------- Weighted data -------------------------
    X1w = X1locentered * L1;     % n x T
    X2w = X2locentered * L2;     % n x T

    % ---------------- Covariances (time metric) -------------
    Sxx_raw = (X1w' * X1w) / n;
    Syy_raw = (X2w' * X2w) / n;
    Sxy_raw = (X1w' * X2w) / n;

    % ---------------- Ridge regularization ------------------
    ridge = ridgeTol;
    Sxx_raw = Sxx_raw + ridge * (trace(Sxx_raw)/T) * eye(T);
    Syy_raw = Syy_raw + ridge * (trace(Syy_raw)/T) * eye(T);

    ResX = stableSigmaSVD(Sxx_raw, svdTol);
    ResY = stableSigmaSVD(Syy_raw, svdTol);
    Sxx = ResX.Sigma;    SxxInv = ResX.SigmaInv;
    Syy = ResY.Sigma;    SyyInv = ResY.SigmaInv;
    Sxy = Sxy_raw;

    % ---------------- Containers & helpers ------------------
    Ak = zeros(T, numOfCCA);  Bk = zeros(T, numOfCCA);  rhoVec = zeros(1, numOfCCA);
    normM = @(v,M) sqrt(max(v' * M * v, eps));

    % ---------------- First pair (alt updates) ---------------
    rng('default');  
    a = randn(T,1); a = a / normM(a,Sxx);
    b = randn(T,1); b = b / normM(b,Syy);
    for it = 1:maxIter
        a_old = a; b_old = b;
        b = SyyInv * (Sxy' * a);  b = b / normM(b,Syy);
        a = SxxInv * (Sxy  * b);  a = a / normM(a,Sxx);
        if norm(a-a_old)+norm(b-b_old) < tol, break; end
    end
    rho = (a' * Sxy * b) / (normM(a,Sxx) * normM(b,Syy));
    Ak(:,1) = a;  Bk(:,1) = b;  rhoVec(1) = rho;

    % ---------------- Higher pairs: deflation  -
    for k = 2:numOfCCA
        Sxy_cur = Sxy;
        for i = 1:k-1
            Sxy_cur = Sxy_cur - (Sxx * Ak(:,i)) * (rhoVec(i)) * (Bk(:,i)' * Syy);
        end
        a = randn(T,1); a = a / normM(a,Sxx);
        b = randn(T,1); b = b / normM(b,Syy);
        for it = 1:maxIter
            a_old = a; b_old = b;
            b = SyyInv * (Sxy_cur' * a);
            if k > 2
                Bprev = Bk(:,1:k-1);
                proj = Bprev' * (Syy * b);
                if norm(proj) > 1e-8 * norm(b), b = b - Bprev * proj; end
            end
            b = b / normM(b,Syy);

            a = SxxInv * (Sxy_cur * b);
            if k > 2
                Aprev = Ak(:,1:k-1);
                proj = Aprev' * (Sxx * a);
                if norm(proj) > 1e-8 * norm(a), a = a - Aprev * proj; end
            end
            a = a / normM(a,Sxx);
            if norm(a-a_old)+norm(b-b_old) < tol, break; end
        end
        rho = (a' * Sxy_cur * b) / (normM(a,Sxx) * normM(b,Syy));
        Ak(:,k) = a;  Bk(:,k) = b;  rhoVec(k) = rho;
    end

    % ---------------- Map directions back to raw time axis ---
    % Scores use X*_w * weights; raw-axis “loadings” use L1/L2 * weights
    Ak_raw = L1 * Ak;    % T x K
    Bk_raw = L2 * Bk;    % T x K


    comp1=X1locentered*Ak_raw;  % n x K
    comp2=X2locentered*Bk_raw;  % n x K


    normX = sqrt(sum(X1locentered.^2, 1))';
    normS = sqrt(sum(comp1.^2, 1));
    corrMat = abs(X1locentered' * comp1) ./ (normX * normS);  % [T x K]
    w_new = sum(corrMat, 2);
    w_new=T*(w_new/sum(w_new));
    

    % ---------------- Pack results ---------------------------
    results.Ak = Ak; results.Bk = Bk;
    results.Ak_raw = Ak_raw; results.Bk_raw = Bk_raw;
    results.rhoVec = rhoVec;

    results.X1_localMean = m1_loc; results.X2_localMean = m2_loc;
    results.X1c = X1locentered; results.X2c = X2locentered;
    results.X1w = X1w; results.X2w = X2w;

    results.comp1=comp1;   results.comp2=comp2; results.Wg=diag(w_new);
    
    results.Sxx = Sxx; results.Syy = Syy; results.Sxy = Sxy;
    results.SxxInv = SxxInv; results.SyyInv = SyyInv;
    results.Sxx_raw = Sxx_raw; results.Syy_raw = Syy_raw; results.Sxy_raw = Sxy_raw;

    results.mixOrder = lower(mixOrder);
    results.L1 = L1; results.L2 = L2;  % keep the transforms for reference

    results.effectiveWeight=sum(L1, 1)';  %col-sum [T x 1];

end



function M = normalize_rows(M)
    rs = sum(M,2); rs(rs==0) = 1; M = M ./ rs;
end


function R = stableSigmaSVD(Sigma, eps_sv)
    Sigma = (Sigma + Sigma')/2;
    [U,S,~] = svd(Sigma);
    s = diag(S);
    s_trunc = s; s_trunc(s <= eps_sv) = 0;
    s_inv = zeros(size(s)); nz = s > eps_sv; s_inv(nz) = 1./s(nz);
    R.Sigma    = U * diag(s_trunc) * U';
    R.SigmaInv = U * diag(s_inv)   * U';
end
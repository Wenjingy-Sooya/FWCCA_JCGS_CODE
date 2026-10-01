function out = EV_from_projection_C(E, G, ridge)
% EV_from_projection_C
%   Project ground-truth maps G onto the span of estimated maps E
%   using least-squares with optional ridge regularization.
%
% Inputs
%   E     : [N x Khat]  estimated spatial maps (columns)
%   G     : [N x K]     ground-truth maps (columns)
%   ridge : scalar ridge parameter (default 1e-8)
%
% Outputs (struct)
%   .C     : [Khat x K] coefficient matrix, Ghat = E * C
%   .Ghat  : [N x K]    projections of G onto span(E)
%   .EV    : [K x 1]    explained variance per GT map = ||ĝ||^2 / ||g||^2
%   .RCOR  : [K x 1]    corr(g, ĝ)
%   .RMSE  : [K x 1]    RMSE per GT map
%   .rankE : numeric    numerical rank of E (after normalization)

    if nargin < 3
        ridge = 1e-8;
    end

    % --- checks ---
    [N, Khat] = size(E);
    [N2, K]   = size(G);
    assert(N == N2, 'E and G must have same number of pixels (rows).');

    % --- column-normalize both E and G ---
    E = E ./ max(vecnorm(E, 2, 1), eps);
    G = G ./ max(vecnorm(G, 2, 1), eps);   % then EV = ||Ghat||^2

    % --- coefficients C = (E'E + rI)^(-1) E' G ---
    Gram = (E' * E) + ridge * eye(Khat);
    C    = Gram \ (E' * G);          % [Khat x K]

    % --- projections ---
    Ghat = E * C;                    % [N x K]

    % --- metrics ---
    gNorm2   = sum(G.^2, 1);                 % 1 × K (||g||²)
    EV       = sum(Ghat.^2, 1) ./ gNorm2;    % 1 × K (||ĝ||² / ||g||²)
    EV       = min(max(EV, 0), 1);           % clip numeric drift
    EV       = EV(:);                        % K × 1

    RCOR = zeros(K, 1);
    RMSE = zeros(K, 1);
    for k = 1:K
        g  = G(:,k);
        gh = Ghat(:,k);
        RCOR(k) = (g' * gh) / (norm(g) * norm(gh) + eps);
        RMSE(k) = norm(g - gh) / sqrt(N);
    end

    % --- numerical rank (for info) ---
    rankE = rank(E);

    % --- pack ---
    out.C     = C;
    out.Ghat  = Ghat;
    out.EV    = EV;
    out.RCOR  = RCOR;
    out.RMSE  = RMSE;
    out.rankE = rankE;
end
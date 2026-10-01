function angles = computeCanonicalAngles(ResultCell, Kcomp)
% computeCanonicalAngles
%
% Computes principal angles between canonical directions obtained under
% different initialization levels.
%
% Input
% -----
% ResultCell : M x 4 cell array
%              columns correspond to
%    
%
% Kcomp      : number of canonical components.
%
% Output
% ------
% angles(m,c,k)
%
% c = 1 : oracle vs noisyOracle
% c = 2 : oracle vs Preliminary CCA
% c = 3 : oracle vs Random matched
%
% k = canonical component
%
% NaN indicates one (or both) runs did not converge.

[M,nInit] = size(ResultCell);

if nInit ~= 4
    error('ResultCell should have four initialization settings.');
end

angles = NaN(M,3,Kcomp);

referenceIdx = 1;
compareIdx = [2 3 4];

for c = 1:3

    idx2 = compareIdx(c);

    for m = 1:M

        R1 = ResultCell{m,referenceIdx};
        R2 = ResultCell{m,idx2};

        %% Require both converged
        if ~(R1.converged && R2.converged)
            continue;
        end

        A1 = R1.Ak_raw;
        A2 = R2.Ak_raw;

        %% Normalize each canonical direction
        A1 = A1 ./ vecnorm(A1);
        A2 = A2 ./ vecnorm(A2);

        %% Match component ordering
        similarity = abs(A1' * A2);

        scoreSame = trace(similarity);

        scoreSwap = 0;
        if Kcomp == 2
            scoreSwap = similarity(1,2) + similarity(2,1);
        end

        if Kcomp == 2 && scoreSwap > scoreSame
            A2 = A2(:,[2 1]);
        end

        %% Compute angles
        for k = 1:Kcomp

            cosineValue = abs(A1(:,k)'*A2(:,k));

            cosineValue = min(max(cosineValue,-1),1);

            angles(m,c,k) = acosd(cosineValue);

        end

    end

end

end
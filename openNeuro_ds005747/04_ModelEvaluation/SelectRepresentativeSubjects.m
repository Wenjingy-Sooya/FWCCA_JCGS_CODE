% =========================================================================
% SelectRepresentativeSubjects.m
%
% Purpose:
%   Select a representative subject for the component-wise improvements
%   observed for the local FWCCA variants in the held-out analysis.
%
%   Based on the group-level component-wise results, the analysis focuses
%   on the components showing the clearest improvements of LG-FWCCA and
%   GL-FWCCA over CCA:
%
%       AnG  : component 4
%       PCgG : component 2
%       PCu  : component 1
%
%   SFG is not included because the local variants do not show a
%   consistent component-wise improvement in this ROI.
%
%   For each selected ROI/component, the LG-FWCCA and GL-FWCCA gains
%   relative to CCA are computed for each subject. The representative
%   subject is defined as the subject whose pair of gains
%
%       (LG-FWCCA - CCA, GL-FWCCA - CCA)
%
%   is closest, in Euclidean distance, to the corresponding group-average
%   pair of gains.
%
% Input:
%   Results/ComponentWiseCorrResults.mat
%
% Output:
%   Results/RepresentativeSubjectResults.mat
%
% =========================================================================
clear;
clc;
close all;

%% ============================================================
% Project paths
%% ============================================================

% Directory containing this script:
%   openNeuro_ds005747/04_ModelEvaluation
scriptDir = fileparts(mfilename('fullpath'));

% Output/results directory
resultsDir = fullfile(scriptDir, 'Results');

%% ============================================================
% Load stored component-wise correlations
%% ============================================================

resultFile = fullfile( ...
    resultsDir, ...
    'ComponentWiseCorrResults.mat');

S = load(resultFile);
ComponentWiseCorrResults = S.ComponentWiseCorrResults;
fprintf('Loaded:\n%s\n', resultFile);


%% ============================================================
% Settings
%% ============================================================

roi_list = {'AnG','PCgG', 'PCu'};
k_list = [4 2 1];
RepresentativeSubjectResults = struct();

%% ============================================================
% Loop over ROIs
%% ============================================================

for roi_idx = 1:numel(roi_list)

    roiName = roi_list{roi_idx};
    k =  k_list(roi_idx);
    R =  ComponentWiseCorrResults.(roiName);
    subjects = R.subjects;
    nSub = numel(subjects);

    %% --------------------------------------------------------
    % Method row indices
    %% --------------------------------------------------------
    methodNames = R.methodNames;
    idxCCA = find(strcmp(methodNames,'CCA'));
    idxLG = find(strcmp(methodNames,'LG-FWCCA'));
    idxGL = find(strcmp(methodNames,'GL-FWCCA'));
    %% --------------------------------------------------------
    % Extract component-wise absolute correlations
    %
    % absCorr:
    %   subject x method x component
    %% --------------------------------------------------------

    absCorr =  R.absCorr;
    rCCA =  squeeze(absCorr(:,idxCCA,k));
    rLG = squeeze(absCorr(:,idxLG,k));
    rGL = squeeze(absCorr(:,idxGL,k));

    %% --------------------------------------------------------
    % Gain relative to CCA
    %% --------------------------------------------------------
    deltaLG =  rLG-rCCA;
    deltaGL =  rGL-rCCA;

    % Group-average gains
    meanDeltaLG = mean(deltaLG,'omitnan');

    meanDeltaGL =  mean(deltaGL,'omitnan');


    %% --------------------------------------------------------
    % Joint distance to group-average gain
    %
    % Each subject:
    %   (deltaLG, deltaGL)
    %
    % Group center:
    %   (meanDeltaLG, meanDeltaGL)
    %% --------------------------------------------------------

    distanceToMean =  sqrt( ...
            (deltaLG-meanDeltaLG).^2 + ...
            (deltaGL-meanDeltaGL).^2);
    %% --------------------------------------------------------
    % Representative subject
    %% --------------------------------------------------------
    [minDistance,repSub] =  min(distanceToMean);
    %% ========================================================
    % Display
    %% ========================================================

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('ROI = %s | component k = %d\n', ...
        roiName,k);
    fprintf('============================================================\n');
    fprintf('Group mean LG gain = %.4f\n', ...
        meanDeltaLG);
    fprintf('Group mean GL gain = %.4f\n\n', ...
        meanDeltaGL);


    for s = 1:nSub

        fprintf( ...
            '%s: CCA=%.3f | LG=%.3f | GL=%.3f | dLG=%+.3f | dGL=%+.3f | distance=%.4f\n', ...
            subjects{s}, ...
            rCCA(s), ...
            rLG(s), ...
            rGL(s), ...
            deltaLG(s), ...
            deltaGL(s), ...
            distanceToMean(s));

    end


    fprintf('\nRepresentative subject = %s\n', ...
        subjects{repSub});

    fprintf('Minimum distance = %.4f\n', ...
        minDistance);


    %% ========================================================
    % Store
    %% ========================================================

    RepresentativeSubjectResults.(roiName).component = k;
    RepresentativeSubjectResults.(roiName).subjectIndex =repSub;
    RepresentativeSubjectResults.(roiName).subject =  subjects{repSub};
    RepresentativeSubjectResults.(roiName).CCA =  rCCA;
    RepresentativeSubjectResults.(roiName).LG = rLG;
    RepresentativeSubjectResults.(roiName).GL = rGL;
    RepresentativeSubjectResults.(roiName).deltaLG = deltaLG;
    RepresentativeSubjectResults.(roiName).deltaGL = deltaGL;
    RepresentativeSubjectResults.(roiName).meanDeltaLG =  meanDeltaLG;
    RepresentativeSubjectResults.(roiName).meanDeltaGL =  meanDeltaGL;
    RepresentativeSubjectResults.(roiName).distanceToMean = distanceToMean;

end

%% ============================================================
% Save representative-subject results
%% ============================================================

saveFile = fullfile( ...
    resultsDir, ...
    'RepresentativeSubjectResults.mat');

save( ...
    saveFile, ...
    'RepresentativeSubjectResults');

fprintf('\nSaved representative-subject results to:\n%s\n', ...
    saveFile);
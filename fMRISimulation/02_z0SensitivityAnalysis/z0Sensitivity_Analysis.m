%% ============================================================
% z0Sensitivity_Analysis for FWCCA
%
% Reviewer comment:
% "Include a simulation study that examines the influence of various
%  methods for obtaining the initial activation indicator z0."
%
% Purpose
% -------
% This simulation investigates the sensitivity of FWCCA to the choice of
% the initial activation indicator z0. Four initialization strategies are
% compared:
%
%   1. Oracle
%   2. Noisy oracle
%   3. Preliminary CCA
%   4. Random matched
%
% Bandwidth specification
% -----------------------
% The local bandwidth h is NOT re-selected within this z0 sensitivity
% experiment. Instead, the fixed bandwidth used at each SNR level is based
% on the preceding bandwidth-selection simulation with 500 Monte Carlo
% replicates.
%
% For each SNR level, we examined the bandwidths selected across the
% 500 simulation replicates and used the bandwidth that was selected most
% frequently (or represented the dominant selection) as the fixed h for
% this sensitivity analysis.
%
% This gives:
%
%   SNR = -5, -4, -3 dB  ->  h = 3.8
%   SNR = -2, ..., 5 dB  ->  h = 3.4
%
% Fixing h in this way isolates the effect of the z0 initialization:
% differences observed in this experiment are therefore attributable to
% the construction of z0 rather than to simultaneous bandwidth tuning.
%
% Simulation design
% -----------------
%   Monte Carlo replicates:  M = 500
%   SNR levels:             -5:5 dB
%   Canonical components:    K = 2
%
% Output
% ------
%   Results/FWCCA_z0Sensitivity_AllSNR.mat
%
% The output file contains:
%   - z0SensitivityResults
%   - summaryTable
%
% The output stores the SNR-specific fixed bandwidth together with the
% corresponding simulation results.


%% ============================================================
clear; clc;

%% ============================================================
% Project directories
%% ============================================================

% Directory containing this script:
%   fMRISimulation/02_z0SensitivityAnalysis/
script_dir = fileparts(mfilename('fullpath'));

% Project root directory:
%   fMRISimulation/
project_root = fileparts(script_dir);

% Shared FWCCA functions:
%   fMRISimulation/Functions/
functions_dir = fullfile(project_root, 'Functions');

if ~exist(functions_dir, 'dir')
    error('Functions directory does not exist:\n%s', functions_dir);
end

% Add the project root so that the MATLAB package
%   fMRISimulation/+z0sensitivity/
% can be accessed as z0sensitivity.*
addpath(project_root);

% Add shared FWCCA functions
addpath(functions_dir);


%% ============================================================
% Output directory
%% ============================================================

% z0 sensitivity-analysis results:
%   fMRISimulation/02_z0SensitivityAnalysis/Results/
output_dir = fullfile(script_dir, 'Results');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end



%% ============================================================
% Basic settings
%% ============================================================
M = 500;
SNR_dB_list = -5:5;
% SNR-dependent fixed bandwidths obtained from the preceding
% 500-replicate bandwidth-selection simulation.
fixed_h_list = [repmat(3.8,1,3) repmat(3.4,1,8)];
grid_size = [20,20];
num_voxels = prod(grid_size);
T = 200;
Kcomp = 2;
tol = 1e-5;
maxIter = 200;
ridge = 1e-5;
svdTol=1e-2;
ridgeTol = 1e-5;
baseSeed = 700000;
maxPlacementTrials = 500;

% z0 initialization methods
z0Methods = {'oracle','noisyOracle','dataDrivenCCA','randomMatched'};
nZ0 = numel(z0Methods);

% Check that each SNR has a corresponding fixed bandwidth
if numel(SNR_dB_list) ~= numel(fixed_h_list)
    error('SNR_dB_list and fixed_h_list must have the same length.');
end

%% ============================================================
% Generate canonical temporal signals
%% ============================================================
boxcar_signal = repmat([ones(1,10),zeros(1,10)],1,T/20);
boxcar_signal = (boxcar_signal-mean(boxcar_signal))/std(boxcar_signal);
quadratic_trend = linspace(-1,1,T).^2;
quadratic_trend = (quadratic_trend-mean(quadratic_trend))/std(quadratic_trend);
%% ============================================================
% Spatial-region settings
%% ============================================================
centerB0 = [9,9];
szB = [6,6];
centerQ0 = [13,13];
szQ = [3,3];
jitRange = [-2,2;-2,2];
I_T = eye(T);
%% ============================================================
% Initialize all-SNR result structure
%% ============================================================
z0SensitivityResults = struct([]);
%% ============================================================
% Main simulation loop
%% ============================================================

for sidx = 1:numel(SNR_dB_list)

    SNR_dB = SNR_dB_list(sidx);
    fixed_h = fixed_h_list(sidx);
    sigma_noise = sqrt(1/10^(SNR_dB/10));
    % Local kernel corresponding to the fixed bandwidth at this SNR
    K_h = makeLocalKernel(T,5,'gaussian',fixed_h);
    fprintf('\n========================================\n');
    fprintf('SNR = %+d dB, fixed h = %.1f\n',SNR_dB,fixed_h);
    fprintf('========================================\n');
    %% --------------------------------------------------------
    % Result containers for current SNR
    %% --------------------------------------------------------
    GFCCAresults_currentSNR = cell(M,nZ0);
    LGFCCAresults_currentSNR = cell(M,nZ0);
    GLFCCAresults_currentSNR = cell(M,nZ0);
    %% ========================================================
    % Monte Carlo loop
    %% ========================================================

    for m = 1:M

        fprintf('SNR %+d dB, run %d / %d\n',SNR_dB,m,M);
        %% ----------------------------------------------------
        % Generate one simulated dataset
        %% ----------------------------------------------------
        % Unique reproducible random seed for each SNR and replicate
        rng(baseSeed+(sidx-1)*M+m,'twister');
        %% ----------------------------------------------------
        % Place boxcar region
        %% ----------------------------------------------------
        dB = draw_jitter(jitRange);
        centerB = centerB0+dB;
        [rowsB,colsB] = block_from_center(centerB,szB,grid_size);
        [RB,CB] = meshgrid(rowsB,colsB);
        boxcar_idx = sub2ind(grid_size,RB(:),CB(:));
        %% ----------------------------------------------------
        % Place quadratic region with required overlap
        %% ----------------------------------------------------
        placed = false;
        placementTrial = 0;

        while ~placed && placementTrial < maxPlacementTrials

            placementTrial = placementTrial+1;
            dQ = draw_jitter(jitRange);
            centerQ = centerQ0+dQ;
            [rowsQ,colsQ] = block_from_center(centerQ,szQ,grid_size);
            [RQ,CQ] = meshgrid(rowsQ,colsQ);
            quadratic_idx = sub2ind(grid_size,RQ(:),CQ(:));
            placed = ~isempty(intersect(boxcar_idx,quadratic_idx));

        end


        if ~placed
            error('Failed to place overlapping regions at SNR index %d, simulation %d.',sidx,m);
        end


        %% ----------------------------------------------------
        % Generate training and test data
        %% ----------------------------------------------------
        Xsignal = zeros(T,num_voxels);
        Xsignal(:,boxcar_idx) = Xsignal(:,boxcar_idx)+boxcar_signal';
        Xsignal(:,quadratic_idx) = Xsignal(:,quadratic_idx)+quadratic_trend';
        Xtrain = Xsignal+sigma_noise*randn(size(Xsignal));
        Xtest = Xsignal+sigma_noise*randn(size(Xsignal));
        Ytrain = Construct4SpatialNeighborsfmri(Xtrain,grid_size);
        XTrain = Xtrain';
        YTrain = Ytrain';
        XTest = Xtest';
        %% ----------------------------------------------------
        % Locally centered data for z0 construction
        %% ----------------------------------------------------
        
        Xc_for_z0 = XTrain-mean(XTrain,1);
        Yc_for_z0 = YTrain-mean(YTrain,1);

        %% ----------------------------------------------------
        % Ground-truth spatial maps
        %% ----------------------------------------------------

        n = size(XTrain,1);
        G = zeros(n,Kcomp);
        G(boxcar_idx,1) = 1;
        G(quadratic_idx,2) = 1;
        Z_true = sum(G,2)>0;

        %% ----------------------------------------------------
        % Construct four z0 indicators
        %% ----------------------------------------------------

        Z0_all = cell(nZ0,1);

        for zidx = 1:nZ0

            methodName = z0Methods{zidx};

            switch methodName

                case 'oracle'

                    Z0 = Z_true;

                case 'noisyOracle'

                    fpRate = 0.20;
                    fnRate = 0.20;
                    Z0 = z0sensitivity.makeNoisyIndicator(Z_true,fpRate,fnRate);

                case 'dataDrivenCCA'

                    activeFrac = mean(Z_true);
                    Z0 = z0sensitivity.makeDataDrivenZ0_CCA(Xc_for_z0,Yc_for_z0,activeFrac);

                case 'randomMatched'

                    Z0 = z0sensitivity.makeRandomMatchedZ0(Z_true);

                otherwise

                    error('Unknown z0 method: %s',methodName);

            end

            Z0_all{zidx} = Z0(:);

        end

        %% ----------------------------------------------------
        % Run FWCCA for each z0 initialization
        %% ----------------------------------------------------

        for zidx = 1:nZ0

            Z0 = Z0_all{zidx};
            % Initial global weights from z0
            w0 = z0sensitivity.makeInitialWeightsFromZ0(Xc_for_z0,Z0,T);
            Wg_0 = diag(w0);
            %% ------------------------------------------------
            % Run G-, LG-, and GL-FWCCA
            %% ------------------------------------------------
            resG_T = Run_FWCCA_Family(XTrain,YTrain,Wg_0,Wg_0,I_T,I_T,Kcomp,maxIter, tol, 'postScale', ridgeTol,  svdTol);
            resLG_T = Run_FWCCA_Family(XTrain,YTrain,Wg_0,Wg_0,K_h,K_h,Kcomp,maxIter,tol,'postScale', ridgeTol,  svdTol);
            resGL_T = Run_FWCCA_Family(XTrain,YTrain,Wg_0,Wg_0,K_h,K_h,Kcomp,maxIter,tol,'preScale',  ridgeTol,   svdTol);
            %% ------------------------------------------------
            % Apply estimated directions to test data
            %% ------------------------------------------------
            XtestG_locentered = XTest-repmat(resG_T.X1_localMean,n,1);
            XtestLG_locentered = XTest-repmat(resLG_T.X1_localMean,n,1);
            XtestGL_locentered = XTest-repmat(resGL_T.X1_localMean,n,1);
            comp1_test_G = XtestG_locentered*resG_T.Ak_raw;
            comp1_test_LG = XtestLG_locentered*resLG_T.Ak_raw;
            comp1_test_GL = XtestGL_locentered*resGL_T.Ak_raw;
            %% ------------------------------------------------
            % EV and RCOR
            %% ------------------------------------------------
            EV_G_train = EV_from_projection_C(resG_T.comp1,G,ridge);
            EV_G_test = EV_from_projection_C(comp1_test_G,G,ridge);
            EV_LG_train = EV_from_projection_C(resLG_T.comp1,G,ridge);
            EV_LG_test = EV_from_projection_C(comp1_test_LG,G,ridge);
            EV_GL_train = EV_from_projection_C(resGL_T.comp1,G,ridge);
            EV_GL_test = EV_from_projection_C(comp1_test_GL,G,ridge);
            %% ------------------------------------------------
            % AUC
            %% ------------------------------------------------
            auc_G_train = computeAUC(sum(abs(resG_T.comp1),2),Z_true);
            auc_G_test = computeAUC(sum(abs(comp1_test_G),2),Z_true);
            auc_LG_train = computeAUC(sum(abs(resLG_T.comp1),2),Z_true);
            auc_LG_test = computeAUC(sum(abs(comp1_test_LG),2),Z_true);
            auc_GL_train = computeAUC(sum(abs(resGL_T.comp1),2),Z_true);
            auc_GL_test = computeAUC(sum(abs(comp1_test_GL),2),Z_true);
            %% ------------------------------------------------
            % Store G-FWCCA results
            %% ------------------------------------------------
            GFCCAresults_currentSNR{m,zidx} = struct( ...
                'InitialZ0',z0Methods{zidx}, ...
                'AUC_Train',auc_G_train, ...
                'AUC_Test',auc_G_test, ...
                'EV_Train',EV_G_train.EV(:)', ...
                'RCOR_Train',EV_G_train.RCOR(:)', ...
                'EV_Test',EV_G_test.EV(:)', ...
                'RCOR_Test',EV_G_test.RCOR(:)', ...
                'Comp1_Train',resG_T.comp1, ...
                'Comp1_Test',comp1_test_G, ...
                'numIter',resG_T.outerIter, ...
                'converged',resG_T.converged, ...
                'finalError',resG_T.finalError, ...
                'Ak_raw',resG_T.Ak_raw, ...
                'Bk_raw',resG_T.Bk_raw, ...
                'G_true',G, ...
                'Z_true',Z_true, ...
                'Z0',Z0, ...
                'AkrawAngleHistory',resG_T.AkrawAngleHistory);
            %% ------------------------------------------------
            % Store LG-FWCCA results
            %% ------------------------------------------------
            LGFCCAresults_currentSNR{m,zidx} = struct( ...
                'InitialZ0',z0Methods{zidx}, ...
                'AUC_Train',auc_LG_train, ...
                'AUC_Test',auc_LG_test, ...
                'EV_Train',EV_LG_train.EV(:)', ...
                'RCOR_Train',EV_LG_train.RCOR(:)', ...
                'EV_Test',EV_LG_test.EV(:)', ...
                'RCOR_Test',EV_LG_test.RCOR(:)', ...
                'Comp1_Train',resLG_T.comp1, ...
                'Comp1_Test',comp1_test_LG, ...
                'numIter',resLG_T.outerIter, ...
                'converged',resLG_T.converged, ...
                'finalError',resLG_T.finalError, ...
                'Ak_raw',resLG_T.Ak_raw, ...
                'Bk_raw',resLG_T.Bk_raw, ...
                'G_true',G, ...
                'Z_true',Z_true, ...
                'Z0',Z0, ...
                'AkrawAngleHistory',resLG_T.AkrawAngleHistory);
            %% ------------------------------------------------
            % Store GL-FWCCA results
            %% ------------------------------------------------
            GLFCCAresults_currentSNR{m,zidx} = struct( ...
                'InitialZ0',z0Methods{zidx}, ...
                'AUC_Train',auc_GL_train, ...
                'AUC_Test',auc_GL_test, ...
                'EV_Train',EV_GL_train.EV(:)', ...
                'RCOR_Train',EV_GL_train.RCOR(:)', ...
                'EV_Test',EV_GL_test.EV(:)', ...
                'RCOR_Test',EV_GL_test.RCOR(:)', ...
                'Comp1_Train',resGL_T.comp1, ...
                'Comp1_Test',comp1_test_GL, ...
                'numIter',resGL_T.outerIter, ...
                'converged',resGL_T.converged, ...
                'finalError',resGL_T.finalError, ...
                'Ak_raw',resGL_T.Ak_raw, ...
                'Bk_raw',resGL_T.Bk_raw, ...
                'G_true',G, ...
                'Z_true',Z_true, ...
                'Z0',Z0, ...
                'AkrawAngleHistory',resGL_T.AkrawAngleHistory);

        end


    end

    %% ========================================================
    % Store current SNR results
    %% ========================================================
    z0SensitivityResults(sidx).SNR_dB = SNR_dB;
    z0SensitivityResults(sidx).fixed_h = fixed_h;
    z0SensitivityResults(sidx).GFCCAresults_currentSNR = GFCCAresults_currentSNR;
    z0SensitivityResults(sidx).LGFCCAresults_currentSNR = LGFCCAresults_currentSNR;
    z0SensitivityResults(sidx).GLFCCAresults_currentSNR = GLFCCAresults_currentSNR;
    fprintf('Completed SNR = %+d dB with h = %.1f.\n',SNR_dB,fixed_h);

end


%% ============================================================
% Summary of SNR-dependent bandwidths
%% ============================================================
summaryTable = table(SNR_dB_list(:),fixed_h_list(:), ...
    'VariableNames',{'SNR_dB','fixed_h'});

disp(summaryTable);

%% ============================================================
% Save results
%% ============================================================

outputMat = fullfile( ...
    output_dir, ...
    'FWCCA_z0Sensitivity_AllSNR.mat');

save(outputMat, ...
    'z0SensitivityResults', ...
    'summaryTable', ...
    '-v7.3');

fprintf('\nSaved all SNR results to:\n%s\n', outputMat);
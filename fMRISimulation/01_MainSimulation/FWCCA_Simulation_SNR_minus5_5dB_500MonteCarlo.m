%% ============================================================
% FWCCA_Simulation_SNR_minus5_5dB_500MonteCarlo.m
%
% This script reproduces the simulation study used to compare
% CCA, G-FWCCA, LG-FWCCA, and GL-FWCCA.
%
% The simulation uses:
%   - a 20-by-20 spatial grid;
%   - 200 time points;
%   - two underlying temporal signals (boxcar and quadratic);
%   - SNR levels from -5 dB to 5 dB;
%   - 500 Monte Carlo replicates per SNR level;
%   - two canonical components;
%   - a range of Gaussian kernel bandwidths for LG-FWCCA
%     and GL-FWCCA.
%
% For each method, the script evaluates:
%   - training and testing AUC;
%   - explained variation (EV);
%   - reconstruction correlation (RCOR);
%   - local spatial variance of the estimated canonical maps.
%
% Independent noisy training and testing datasets are generated
% for each Monte Carlo replicate. Global temporal weights are
% initialized using a noisy version of the true active-voxel labels.
%
% Output:
%   The simulation results are saved as:
%   Results/FWCCA_results_SNR_minus5_5_Gaussian_5neighborhood.mat
%
%% ============================================================

clear; clc; close all;

%% ============================================================
% Project directories
%% ============================================================

% Directory containing this script:
%   fMRISimulation/01_MainSimulation/
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

% Add shared FWCCA functions to the MATLAB path
addpath(functions_dir);

%% ============================================================
% Output directories
%% ============================================================

% Numerical simulation results:
%   fMRISimulation/01_MainSimulation/Results/
results_dir = fullfile(script_dir, 'Results');

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Simulation figures:
%   fMRISimulation/01_MainSimulation/Figures/
figures_dir = fullfile(script_dir, 'Figures');

if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
%% ============================================================
% Configuration
%% ============================================================
% Simulations per SNR 
M =500;
% SNR levels (variance-based, dB)
SNR_dB_list = -5:1:5;
% Grid and time
grid_size = [20, 20];
num_voxels = prod(grid_size);
num_timepoints = 200;
Kcomp = 2;                                      % extract 2 comps

nSNR = numel(SNR_dB_list);
nMeth = 2;
h_list = [0.5 1.0 1.8 2.8 3.4 3.8 4.2 4.6 5];   % bandwidths
n_h = numel(h_list);
% Fix dimension ordering: (nMeth, nSNR, M)
LocalWeightVariantresults = cell(nMeth, nSNR, M, n_h);  % NOTE: meth × snr × run

for midx = 1:nMeth      % method index (1=LG, 2=GL)
  for sidx = 1:nSNR     % SNR index
     for m = 1:M        % Monte Carlo repetition
         for h_idx=1:n_h
            S = struct(...
            'LocalVals',    zeros(1, Kcomp), ...
            'AvgLocalVals', zeros(1, 1), ...
            'AUCs_Train',   zeros(1, 1), ...
            'AUCs_Test',    zeros(1, 1), ...
            'EV_Train',     zeros(1, Kcomp), ...
            'RCOR_Train',   zeros(1, Kcomp), ...
            'EV_Test',      zeros(1, Kcomp), ...
            'RCOR_Test',    zeros(1, Kcomp), ...
            'Comp1_Test',   zeros(num_voxels, Kcomp) ...
           );
            LocalWeightVariantresults{midx, sidx, m, h_idx} = S;
         end
     end
  end
end

% Define structure to hold results for CCA/GFCCA 1:CCA; 2:GFWCCA
results = cell(nMeth,nSNR, M); 
for midx = 1:nMeth   
  for sidx = 1:nSNR
     for  m = 1:M
        results{midx,sidx,m} = struct( ...
        'LocalVals',    zeros(1, Kcomp), ...
        'AvgLocalVals', zeros(1, 1),...
        'AUCs_Train',   zeros(1, 1), ...
        'AUCs_Test',    zeros(1, 1), ...
        'EV_Train',     zeros(1, Kcomp), ...     % [EV_boxcar, EV_quadratic] per sim
        'RCOR_Train',   zeros(1, Kcomp), ...     % [RC_boxcar, RC_quadratic] per sim
        'EV_Test',      zeros(1, Kcomp), ...     % [EV_boxcar, EV_quadratic] per sim
        'RCOR_Test',    zeros(1, Kcomp), ...     % [RC_boxcar, RC_quadratic] per sim
        'Comp1_Test',   cell(1,1)...
       );
     end
  end
end



%% Generate canonical signals
boxcar_signal = repmat([ones(1,10), zeros(1,10)], 1, num_timepoints/20);
boxcar_signal = (boxcar_signal - mean(boxcar_signal)) / std(boxcar_signal);

quadratic_trend = (linspace(-1,1,num_timepoints)).^2;
quadratic_trend = (quadratic_trend - mean(quadratic_trend)) / std(quadratic_trend);

centerB0 = [9, 9];    szB = [6, 6];
centerQ0 = [13, 13];  szQ = [3, 3];
jitRange = [-2 2; -2 2];

svdTol   = 1e-1;
ridgeTol = 1e-4;

for sidx = 1:nSNR

  SNR_dB = SNR_dB_list(sidx);
  sigma_noise = sqrt(1 / (10^(SNR_dB/10)));

  for m = 1:M
    rng(100000 * sidx + m);
    % ----- Place Boxcar region (mean + jitter) -----
      dB = draw_jitter(jitRange);
      centerB = centerB0 + dB;
      [rowsB, colsB] = block_from_center(centerB, szB, grid_size);
      [RB, CB] = meshgrid(rowsB, colsB);
      boxcar_idx = sub2ind(grid_size, RB(:), CB(:));

   % ----- Place Quadratic region with overlap policy  -----
    placed = false;
    
    while ~placed
      dQ = draw_jitter(jitRange);
      centerQ = centerQ0 + dQ;
      [rowsQ, colsQ] = block_from_center(centerQ, szQ, grid_size);
      [RQ, CQ] = meshgrid(rowsQ, colsQ);
      quadratic_idx = sub2ind(grid_size, RQ(:), CQ(:));
      
      if ~isempty(intersect(boxcar_idx, quadratic_idx))
        placed = true;
      end

    end
    % Create data
    Xtrain = zeros(num_timepoints, num_voxels);
    Xtrain(:, boxcar_idx) = Xtrain(:, boxcar_idx) + boxcar_signal';
    Xtrain(:, quadratic_idx) = Xtrain(:, quadratic_idx) + quadratic_trend';
    Xtest = Xtrain;

    Xtrain = Xtrain + sigma_noise * randn(size(Xtrain));
    Xtest  = Xtest  + sigma_noise * randn(size(Xtest));

    Ytrain = Construct4SpatialNeighborsfmri(Xtrain, grid_size);
   
    XTrain = Xtrain'; YTrain = Ytrain'; XTest = Xtest'; 

    % Label vector Z
    n = size(XTrain, 1);
    G = zeros(n, 2);
    G(boxcar_idx, 1) = 1;
    G(quadratic_idx, 2) = 1;
    Z_true = sum(G, 2) > 0;

    % Noisy label init
    Z0 = Z_true;
    idx_fp = randsample(find(~Z0), round(0.2 * sum(~Z0)));
    idx_fn = randsample(find(Z0), round(0.2 * sum(Z0)));
    Z0(idx_fp) = 1; Z0(idx_fn) = 0;

    % Compute global weights
    Z0c = Z0 - mean(Z0);
    Xc = XTrain - mean(XTrain,1);
    corrs = abs((Xc' * Z0c) ./ sqrt(sum(Xc.^2,1)' * sum(Z0c.^2)));
    corrs(isnan(corrs)) = 0;
    w0 = num_timepoints * corrs / sum(corrs);
    Wg_0 = diag(w0);



     %% CCA######################################
      midx=1;
      resCCA_T = Run_FWCCA_Family(XTrain, YTrain, eye(num_timepoints), eye(num_timepoints), eye(num_timepoints), eye(num_timepoints), Kcomp, 100, 1e-5, 'postScale', ridgeTol, svdTol);
      XtestCCA_locentered = XTest- repmat(resCCA_T.X1_localMean, n, 1);
      resCCA_E.comp1=XtestCCA_locentered*resCCA_T.Ak_raw;
      results{midx,sidx,m}.Comp1_Test=resCCA_E.comp1;

      EV_CCAT = EV_from_projection_C(resCCA_T.comp1, G, 1e-5);
      EV_CCAE = EV_from_projection_C(resCCA_E.comp1, G, 1e-5);

      results{midx,sidx,m}.EV_Train=EV_CCAT.EV;
      results{midx,sidx,m}.EV_Test=EV_CCAE.EV;
      
      results{midx,sidx,m}.RCOR_Train=EV_CCAT.RCOR;
      results{midx,sidx,m}.RCOR_Test=EV_CCAE.RCOR;

      combinedCCA_T_map = sum(abs(resCCA_T.comp1), 2);  % Combine abs(Comp1) + abs(Comp2)
      results{midx,sidx,m}.AUCs_Train = computeAUC(combinedCCA_T_map, Z_true);
       
      combinedCCA_E_map = sum(abs(resCCA_E.comp1), 2);  % Combine abs(Comp1) + abs(Comp2)
      results{midx,sidx,m}.AUCs_Test = computeAUC(combinedCCA_E_map, Z_true);

      %% G-FWCCA######################################

      midx=2;
      
      resG_T = Run_FWCCA_Family(XTrain, YTrain, Wg_0, Wg_0, eye(num_timepoints), eye(num_timepoints), Kcomp, 100, 1e-5, 'postScale', ridgeTol, svdTol);
    
      XtestG_locentered = XTest- repmat(resG_T.X1_localMean, n, 1);
      resG_E.comp1=XtestG_locentered*resG_T.Ak_raw;
      results{midx,sidx,m}.Comp1_Test=resG_E.comp1;

      EV_GT = EV_from_projection_C(resG_T.comp1, G, 1e-5);
      EV_GE = EV_from_projection_C(resG_E.comp1, G, 1e-5);

      results{midx,sidx,m}.EV_Train=EV_GT.EV;
      results{midx,sidx,m}.EV_Test=EV_GE.EV;
      
      results{midx,sidx,m}.RCOR_Train=EV_GT.RCOR;
      results{midx,sidx,m}.RCOR_Test=EV_GE.RCOR;

      combinedG_T_map = sum(abs(resG_T.comp1), 2);  % Combine abs(Comp1) + abs(Comp2)
      results{midx,sidx,m}.AUCs_Train = computeAUC(combinedG_T_map, Z_true);
       
      combinedG_E_map = sum(abs(resG_E.comp1), 2);  % Combine abs(Comp1) + abs(Comp2)
      results{midx,sidx,m}.AUCs_Test= computeAUC(combinedG_E_map, Z_true);
      
      %calculate localvar for each canonical maps.

      for k = 1:Kcomp
        map_CCA = reshape(resCCA_T.comp1(:,k), grid_size); 
        map_G = reshape(resG_T.comp1(:,k), grid_size);        
        results{1,sidx,m}.LocalVals(k) = computeLocalVariance3x3(map_CCA);
        results{2,sidx,m}.LocalVals(k) = computeLocalVariance3x3(map_G);     
      end


      for midx = 1:nMeth        
           results{midx,sidx,m}.AvgLocalVals=mean(results{midx,sidx,m}.LocalVals(:));
      end
      

    for h_idx = 1:n_h
      
      h = h_list(h_idx);
      K = makeLocalKernel(num_timepoints, 5, 'gaussian', h);
      fprintf('SNR =%.2f, Run %d, bandwidth h = %.2f for local feature weight variants\n',SNR_dB, m, h);
      
      %% LG-FWCCA ######################################
      midx=1;
      resLG_T = Run_FWCCA_Family(XTrain, YTrain, Wg_0, Wg_0, K, K, Kcomp, 100, 1e-5, 'postScale', ridgeTol, svdTol);
      XtestLG_locentered = XTest-repmat(resLG_T.X1_localMean, n, 1);
      
      resLG_E.comp1=XtestLG_locentered*resLG_T.Ak_raw;
      
      LocalWeightVariantresults{midx, sidx, m, h_idx}.Comp1_Test =  resLG_E.comp1;

      EV_LGT = EV_from_projection_C(resLG_T.comp1, G, 1e-5);
      EV_LGE = EV_from_projection_C(resLG_E.comp1, G, 1e-5);

      LocalWeightVariantresults{midx, sidx, m, h_idx}.EV_Train=EV_LGT.EV;
      LocalWeightVariantresults{midx, sidx, m, h_idx}.EV_Test=EV_LGE.EV;
      
      LocalWeightVariantresults{midx, sidx, m, h_idx}.RCOR_Train=EV_LGT.RCOR;
      LocalWeightVariantresults{midx, sidx, m, h_idx}.RCOR_Test=EV_LGE.RCOR;

      combinedLG_T_map = sum(abs(resLG_T.comp1), 2);  % Combine abs(Comp1) + abs(Comp2)
      LocalWeightVariantresults{midx, sidx, m, h_idx}.AUCs_Train=computeAUC(combinedLG_T_map, Z_true);
       
      combinedLG_E_map = sum(abs(resLG_E.comp1), 2);  % Combine abs(Comp1) + abs(Comp2)
      LocalWeightVariantresults{midx, sidx, m, h_idx}.AUCs_Test=computeAUC(combinedLG_E_map, Z_true);

 
      %% GL-FWCCA ######################################
      midx=2;
      resGL_T = Run_FWCCA_Family(XTrain, YTrain, Wg_0, Wg_0, K, K, Kcomp, 100, 1e-5, 'preScale', ridgeTol,  svdTol);
      XtestGL_locentered = XTest- repmat(resGL_T.X1_localMean, n, 1);

      resGL_E.comp1=XtestGL_locentered*resGL_T.Ak_raw;
      LocalWeightVariantresults{midx,sidx,m ,h_idx}.Comp1_Test=resGL_E.comp1;

      EV_GLT = EV_from_projection_C(resGL_T.comp1, G, 1e-5);
      EV_GLE = EV_from_projection_C(resGL_E.comp1, G, 1e-5);


      LocalWeightVariantresults{midx,sidx,m, h_idx}.EV_Train=EV_GLT.EV;
      LocalWeightVariantresults{midx,sidx,m, h_idx}.EV_Test=EV_GLE.EV;
      
      LocalWeightVariantresults{midx,sidx,m, h_idx}.RCOR_Train=EV_GLT.RCOR;
      LocalWeightVariantresults{midx,sidx,m, h_idx}.RCOR_Test=EV_GLE.RCOR;

      combinedGL_T_map = sum(abs(resGL_T.comp1), 2);  % Combine abs(Comp1) + abs(Comp2)
      LocalWeightVariantresults{midx,sidx,m, h_idx}.AUCs_Train = computeAUC(combinedGL_T_map, Z_true);
       
      combinedGL_E_map = sum(abs(resGL_E.comp1), 2);  % Combine abs(Comp1) + abs(Comp2)
      LocalWeightVariantresults{midx,sidx,m, h_idx}.AUCs_Test= computeAUC(combinedGL_E_map, Z_true);
      
      
      %calculate localvar for each canonical maps.

      for k = 1:Kcomp
        map_LG = reshape(resLG_T.comp1(:,k), grid_size); 
        map_GL = reshape(resGL_T.comp1(:,k), grid_size);     
        LocalWeightVariantresults{1,sidx,m,h_idx}.LocalVals(k) = computeLocalVariance3x3(map_LG);
        LocalWeightVariantresults{2,sidx,m,h_idx}.LocalVals(k) = computeLocalVariance3x3(map_GL);
        
      end


      for midx = 1:nMeth 
           LocalWeightVariantresults{midx,sidx,m, h_idx}.AvgLocalVals=mean(LocalWeightVariantresults{midx,sidx,m, h_idx}.LocalVals);
      end
 

     %end for h
    end

    %end for m
  end

 
end

%% ============================================================
% Save simulation results
%% ============================================================

filename = fullfile( ...
    results_dir, ...
    'FWCCA_results_SNR_minus5_5_Gaussian_5neighborhood.mat');

% Save all workspace variables with large-variable support
save(filename, '-v7.3');

fprintf('\n============================================================\n');
fprintf('Simulation completed successfully.\n');
fprintf('Results saved to:\n%s\n', filename);
fprintf('============================================================\n');
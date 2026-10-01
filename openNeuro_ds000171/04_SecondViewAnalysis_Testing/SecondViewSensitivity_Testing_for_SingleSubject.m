%%
% Held-out testing analysis for the FWCCA second-view sensitivity study
% using the OpenNeuro Emotional Music dataset (ds000171).
%
% This script evaluates CCA, G-FWCCA, LG-FWCCA, and GL-FWCCA on the
% independent testrun for one subject. All model fitting, second-view
% construction, and hyperparameter selection are based exclusively on the
% trainrun results generated in 03_SecondViewSensitivityAnalysis.
%
% INPUTS
% ------
% 1. Training inputs:
%    03_SecondViewSensitivityAnalysis/SubjectResults/<group>/<subject>/
%        <subject>_SecondViewInputs.mat
%
%    Contains the training masks, X1/X2 matrices, global feature weights,
%    second-view thresholds, and analysis settings.
%
% 2. CCA / G-FWCCA training results:
%        <subject>_CCA_GFWCCA_TrainingResults.mat
%
%    Contains the CCA and G-FWCCA models fitted using the trainrun.
%
% 3. LG / GL-FWCCA cross-validation results:
%        <subject>_LG_GLFWCCA_CVResults.mat
%
%    Contains the K-specific CV-selected local-kernel hyperparameters for
%    LG-FWCCA and GL-FWCCA.
%
% 4. Held-out testing fMRI:
%    Data/<group>/<subject>/func/
%        swar<subject>_task-music_testrun_bold.nii
%
% 5. Subject-specific SPM model:
%    Data/<group>/<subject>/first_level_analysis/SPM.mat
%
% 6. Functional-system mapping:
%    Atlas_FunctionMap_Multilabel.xlsx
%
% OUTPUTS
% -------
% Subject-level held-out functional-system counts are saved under:
%
%    04_SecondViewAnalysis_Testing/Results/<group>/<subject>/
%        K1/
%            positive_MethodComparison.xlsx
%            negative_MethodComparison.xlsx
%        ...
%        K<Kcomp_max>/
%            positive_MethodComparison.xlsx
%            negative_MethodComparison.xlsx
%
% Each workbook contains one worksheet for each second-view threshold
% (p90, p85, p80, p75, p70, and p65). Each worksheet contains:
%
%    System | CCA | G_FWCCA | LG_FWCCA | GL_FWCCA
%
% where the method columns contain cumulative functional-system counts over
% canonical components 1:K.
%
% PROCEDURE
% ---------
% - The training-derived brain mask is applied to the held-out testrun.
% - CCA and G-FWCCA canonical directions fitted on the trainrun are applied
%   directly to the held-out testing data.
% - For LG-FWCCA and GL-FWCCA, the K-specific hyperparameters selected by
%   cross-validation are used to refit the models on the trainrun.
% - The resulting canonical directions are then applied to the testrun.
% - Canonical maps are mapped to Neuromorphometrics functional systems.
% - Functional-system counts are accumulated across components 1:K.
%
% IMPORTANT
% ---------
% The held-out testrun is used only for independent evaluation. It is not
% used for second-view construction, hyperparameter selection, or model
% fitting.
% =========================================================================

%% ============================================================
% Subject settings
% ============================================================

subj = 'sub-mdd02';
groupName = 'MDD';
runName = 'testrun';


%% ============================================================
% Project paths
% ============================================================

testingAnalysisRoot = fileparts(mfilename('fullpath'));
datasetRoot = fileparts(testingAnalysisRoot);
projectRoot = fileparts(datasetRoot);

spm_dir = fullfile(projectRoot, 'spm_25');
functions_dir = fullfile(datasetRoot, 'Functions');

if ~exist(spm_dir, 'dir')
    error('SPM directory does not exist:\n%s', spm_dir);
end

if ~exist(functions_dir, 'dir')
    error('Functions directory does not exist:\n%s', functions_dir);
end


%% ============================================================
% SPM and function setup
% ============================================================

addpath(spm_dir);
addpath(functions_dir);

spm('defaults', 'fmri');
spm_jobman('initcfg');

atlas = spm_atlas('load', 'Neuromorphometrics');


%% ============================================================
% Training / CV result directory
% ============================================================

trainingResultsRoot = fullfile(datasetRoot, '03_SecondViewSensitivityAnalysis', ...
    'SubjectResults', groupName, subj);


%% ============================================================
% Input files
% ============================================================

secondViewInputFile = fullfile(trainingResultsRoot, ...
    sprintf('%s_SecondViewInputs.mat', subj));

globalTrainingFile = fullfile(trainingResultsRoot, ...
    sprintf('%s_CCA_GFWCCA_TrainingResults.mat', subj));

localCVFile = fullfile(trainingResultsRoot, ...
    sprintf('%s_LG_GLFWCCA_CVResults.mat', subj));


%% ============================================================
% Check input files
% ============================================================

if ~exist(secondViewInputFile, 'file')
    error('Cannot find SecondView input file:\n%s', secondViewInputFile);
end

if ~exist(globalTrainingFile, 'file')
    error('Cannot find CCA/G-FWCCA training result:\n%s', globalTrainingFile);
end

if ~exist(localCVFile, 'file')
    error('Cannot find LG/GL-FWCCA CV result:\n%s', localCVFile);
end


%% ============================================================
% Load training inputs and results
% ============================================================

fprintf('\nLoading SecondView inputs:\n%s\n', secondViewInputFile);
S = load(secondViewInputFile, 'AllTrainInputs');
AllTrainInputs = S.AllTrainInputs;
clear S;

fprintf('\nLoading CCA / G-FWCCA training results:\n%s\n', globalTrainingFile);
S = load(globalTrainingFile, 'CCA_GFWCCA_TrainingResults');
CCA_GFWCCA_TrainingResults = S.CCA_GFWCCA_TrainingResults;
clear S;

fprintf('\nLoading LG / GL-FWCCA CV results:\n%s\n', localCVFile);
S = load(localCVFile, 'LG_GLFWCCA_TrainingResults');
LG_GLFWCCA_TrainingResults = S.LG_GLFWCCA_TrainingResults;
clear S;


%% ============================================================
% Analysis settings recovered from training
% ============================================================

Kcomp_max = AllTrainInputs{1}.Kcomp_max;
svdTol = AllTrainInputs{1}.svdTol;
nViews = numel(AllTrainInputs);
ridgeTol = 1e-5;

fprintf('\n');
fprintf('============================================================\n');
fprintf('Testing analysis settings\n');
fprintf('============================================================\n');
fprintf('Subject   : %s\n', subj);
fprintf('Group     : %s\n', groupName);
fprintf('Run       : %s\n', runName);
fprintf('Kcomp_max : %d\n', Kcomp_max);
fprintf('svdTol    : %.6g\n', svdTol);
fprintf('Views     : %d\n', nViews);
fprintf('============================================================\n');


%% ============================================================
% Functional systems
% ============================================================

systems_ordered = {'EmotionControl', 'Visual', 'Auditory', 'Motor', ...
    'LanguageCognition', 'Memory', 'Attention'};

nSystems = numel(systems_ordered);


%% ============================================================
% Load functional-system mapping
% ============================================================

funcmap_file = fullfile(projectRoot, 'Atlas_FunctionMap_Multilabel.xlsx');

if ~exist(funcmap_file, 'file')
    error('Cannot find functional-system mapping file:\n%s', funcmap_file);
end

funcmap = readtable(funcmap_file);
label_names = funcmap.LabelName;
label_systems = funcmap.FunctionalSystem;


%% ============================================================
% Subject paths
% ============================================================

root = fullfile(datasetRoot, 'data', groupName);
subj_dir = fullfile(root, subj);
func_dir = fullfile(subj_dir, 'func');
spm_path = fullfile(subj_dir, 'first_level_analysis', 'SPM.mat');

if ~exist(subj_dir, 'dir')
    error('Subject directory does not exist:\n%s', subj_dir);
end

if ~exist(func_dir, 'dir')
    error('Functional directory does not exist:\n%s', func_dir);
end

if ~exist(spm_path, 'file')
    error('SPM.mat does not exist:\n%s', spm_path);
end


%% ============================================================
% Find held-out TEST functional image
% ============================================================

functional_pattern = ['^swar' subj '_task-music_' runName '_bold\.nii$'];
files = spm_select('FPList', func_dir, functional_pattern);

if isempty(files)
    error('Cannot find testing fMRI for %s.', subj);
end

if size(files,1) > 1
    error('More than one testing functional image was found for %s.', subj);
end

nii_path = strtrim(files(1,:));

fprintf('\nTesting functional image:\n%s\n', nii_path);


%% ============================================================
% Load SPM model and check TR
% ============================================================

load(spm_path, 'SPM');

TR_test = SPM.xY.RT;
TR_train = AllTrainInputs{1}.TR;

fprintf('\nTesting TR = %.3f sec\n', TR_test);

if abs(TR_train - TR_test) > 1e-10
    error(['Training and testing TR values do not match.\n' ...
        'TR_train = %.6f\nTR_test  = %.6f'], TR_train, TR_test);
end

TR = TR_train;


%% ============================================================
% Drift removal
% ============================================================

fprintf('\nRemoving drift from testing data...\n');
Y = get_driftremoved_fMRI(spm_path, nii_path);


%% ============================================================
% Temporal prewhitening
% ============================================================

fprintf('\nApplying temporal prewhitening...\n');

V = full(SPM.xVi.V);

if size(V,1) ~= size(Y,2)
    error(['Temporal covariance dimension mismatch.\n' ...
        'size(V,1) = %d\nsize(Y,2) = %d'], size(V,1), size(Y,2));
end

Vw = sqrtm(inv(V));
Y = Y * Vw;

if ~isreal(Y)
    warning('Complex values produced during whitening. Taking real part.');
    Y = real(Y);
end

T_test = size(Y,2);

fprintf('Testing time points = %d\n', T_test);


%% ============================================================
% Check training / testing temporal dimensions
% ============================================================

T_train = size(AllTrainInputs{1}.X1Mask, 2);

if T_train ~= T_test
    error(['Training and testing time dimensions do not match.\n' ...
        'T_train = %d\nT_test  = %d'], T_train, T_test);
end

fprintf('Training/testing temporal dimensions match: T = %d\n', T_train);


%% ============================================================
% NIfTI header
% ============================================================

Vnii = spm_vol(nii_path);


%% ============================================================
% Output directory
% ============================================================

outputRoot = fullfile(testingAnalysisRoot, 'Results', groupName, subj);

if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

fprintf('\nTesting results will be saved to:\n%s\n', outputRoot);


%% ============================================================
% Create K-specific output directories
% ============================================================

for id_k = 1:Kcomp_max
    K_outputRoot = fullfile(outputRoot, sprintf('K%d', id_k));

    if ~exist(K_outputRoot, 'dir')
        mkdir(K_outputRoot);
    end
end


%% ============================================================
% Basic consistency checks
% ============================================================

if numel(CCA_GFWCCA_TrainingResults) ~= nViews
    error(['Number of views differs between AllTrainInputs and ' ...
        'CCA_GFWCCA_TrainingResults.']);
end

if numel(LG_GLFWCCA_TrainingResults) ~= nViews
    error(['Number of views differs between AllTrainInputs and ' ...
        'LG_GLFWCCA_TrainingResults.']);
end

%% ============================================================
% Loop over positive / negative views
% ============================================================

for viewIdx = 1:nViews

    %% --------------------------------------------------------
    % Current training structures
    %% --------------------------------------------------------

    trainInput = AllTrainInputs{viewIdx};
    globalTrainingOutput = CCA_GFWCCA_TrainingResults{viewIdx};
    localCVOutput = LG_GLFWCCA_TrainingResults(viewIdx);

    viewType = trainInput.viewType;

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('Subject = %s | View = %s\n', subj, viewType);
    fprintf('============================================================\n');


    %% ========================================================
    % Check view correspondence
    % ========================================================

    if isfield(globalTrainingOutput, 'viewType')
        if ~strcmp(globalTrainingOutput.viewType, viewType)
            error(['View mismatch between SecondViewInputs and ' ...
                'CCA/G-FWCCA training results for view %d.'], viewIdx);
        end
    end

    if isfield(localCVOutput, 'viewType')
        if ~strcmp(localCVOutput.viewType, viewType)
            error(['View mismatch between SecondViewInputs and ' ...
                'LG/GL-FWCCA CV results for view %d.'], viewIdx);
        end
    end


    %% ========================================================
    % Training information common to all thresholds
    % ========================================================

    X1Train = trainInput.X1Mask;
    idx_mask = trainInput.idx_mask;
    final_mask = trainInput.final_mask;
    Wg_0 = trainInput.Wg_0;

    vol_size = size(final_mask);
    T_train = size(X1Train, 2);


    %% --------------------------------------------------------
    % Check train/test temporal dimensions
    %% --------------------------------------------------------

    if T_train ~= T_test
        error(['Training and testing time dimensions do not match.\n' ...
            'T_train = %d\nT_test  = %d'], T_train, T_test);
    end


    %% --------------------------------------------------------
    % Extract TEST voxels using training mask
    %% --------------------------------------------------------

    XMask_test = Y(idx_mask, :);

    if size(XMask_test, 1) ~= size(X1Train, 1)
        error(['Training/testing mask voxel counts do not match.\n' ...
            'Training voxels = %d\nTesting voxels  = %d'], ...
            size(X1Train, 1), size(XMask_test, 1));
    end


    %% ========================================================
    % Threshold results
    % ========================================================

    if ~isfield(localCVOutput, 'thresholdResults')
        error('thresholdResults is missing from LG_GLFWCCA_TrainingResults(%d).', viewIdx);
    end

    if ~isfield(globalTrainingOutput, 'results')
        error('results is missing from CCA_GFWCCA_TrainingResults{%d}.', viewIdx);
    end

    LocalFWCCAoutput = localCVOutput.thresholdResults;
    Globaloutput = globalTrainingOutput.results;

    nThresholds = numel(trainInput.secondViewThresholds);


    %% --------------------------------------------------------
    % Check threshold counts
    %% --------------------------------------------------------

    if numel(LocalFWCCAoutput) ~= nThresholds
        error(['Number of thresholds differs between SecondViewInputs and ' ...
            'LG/GL-FWCCA CV results for view %d.'], viewIdx);
    end

    if numel(Globaloutput) ~= nThresholds
        error(['Number of thresholds differs between SecondViewInputs and ' ...
            'CCA/G-FWCCA training results for view %d.'], viewIdx);
    end


    %% ========================================================
    % Delete old workbooks for current view ONCE
    %
    % Different thresholds are written as different worksheets,
    % so the workbook must not be deleted inside the threshold loop.
    % ========================================================

    for id_k = 1:Kcomp_max

        K_outputRoot = fullfile(outputRoot, sprintf('K%d', id_k));
        comparisonWorkbook = fullfile(K_outputRoot, sprintf('%s_MethodComparison.xlsx', viewType));

        if exist(comparisonWorkbook, 'file')
            delete(comparisonWorkbook);
        end

    end


    %% ========================================================
    % Loop over second-view thresholds
    % ========================================================

    for id_p = 1:nThresholds

        %% ----------------------------------------------------
        % Threshold information
        %% ----------------------------------------------------

        percentileValue = LocalFWCCAoutput(id_p).percentile;
        topProportion = LocalFWCCAoutput(id_p).topProportion;


        %% ----------------------------------------------------
        % Current second view
        %
        % X2Mask is stored in SecondViewInputs rather than the
        % cleaned CCA/G-FWCCA training-result file.
        %% ----------------------------------------------------

        X2Train = trainInput.secondViewThresholds(id_p).X2Mask;


        %% ----------------------------------------------------
        % Basic X1 / X2 consistency
        %% ----------------------------------------------------

        if ~isequal(size(X1Train), size(X2Train))
            error(['X1Train and X2Train dimensions do not match.\n' ...
                'View = %s, threshold = %d'], viewType, percentileValue);
        end


        %% ----------------------------------------------------
        % Number of top voxels used for atlas accounting
        %% ----------------------------------------------------

        nVoxels = size(X2Train, 1);
        nTopVoxels = floor(0.0005 * nVoxels);

        fprintf('\n');
        fprintf('------------------------------------------------------------\n');
        fprintf('View = %s | p%d | top %.0f%%\n', viewType, percentileValue, topProportion);
        fprintf('Number of training voxels = %d\n', nVoxels);
        fprintf('nTopVoxels = %d\n', nTopVoxels);
        fprintf('------------------------------------------------------------\n');


        %% ====================================================
        % GLOBAL METHODS
        %
        % CCA and G-FWCCA were already fitted during training.
        % No refitting is performed here.
        % ====================================================

        resCCA_T = Globaloutput(id_p).resCCA_T;
        resG_FWCCA_T = Globaloutput(id_p).resG_FWCCA_T;


        %% ----------------------------------------------------
        % Check number of available components
        %% ----------------------------------------------------

        if size(resCCA_T.Ak_raw, 2) < Kcomp_max
            error('resCCA_T contains fewer than Kcomp_max canonical components.');
        end

        if size(resG_FWCCA_T.Ak_raw, 2) < Kcomp_max
            error('resG_FWCCA_T contains fewer than Kcomp_max canonical components.');
        end


        %% ----------------------------------------------------
        % Project TEST data for CCA and G-FWCCA
        %
        % Output size: Nvox x Kcomp_max
        %% ----------------------------------------------------

        CCAmaps_test = (XMask_test - resCCA_T.X1_localMean) * ...
            resCCA_T.Ak_raw(:, 1:Kcomp_max);

        G_FWCCAmaps_test = (XMask_test - resG_FWCCA_T.X1_localMean) * ...
            resG_FWCCA_T.Ak_raw(:, 1:Kcomp_max);


        %% ====================================================
        % GLOBAL METHOD component-level functional-system counts
        %
        % Compute once for all K.
        % ====================================================

        CCA_component_count = zeros(nSystems, Kcomp_max);
        G_FWCCA_component_count = zeros(nSystems, Kcomp_max);

        for comp = 1:Kcomp_max

            %% ------------------------------------------------
            % CCA
            %% ------------------------------------------------

            canonical_map = CCAmaps_test(:, comp);

            [~, system_table] = CanonicalMap_to_FunctionalSystem( ...
                canonical_map, idx_mask, vol_size, Vnii, atlas, ...
                label_names, label_systems, systems_ordered, nTopVoxels, comp);

            CCA_component_count(:, comp) = system_table.Count;


            %% ------------------------------------------------
            % G-FWCCA
            %% ------------------------------------------------

            canonical_map = G_FWCCAmaps_test(:, comp);

            [~, system_table] = CanonicalMap_to_FunctionalSystem( ...
                canonical_map, idx_mask, vol_size, Vnii, atlas, ...
                label_names, label_systems, systems_ordered, nTopVoxels, comp);

            G_FWCCA_component_count(:, comp) = system_table.Count;

        end


        %% ----------------------------------------------------
        % Cumulative global counts
        %
        % Column K = Comp1 + ... + CompK
        %% ----------------------------------------------------

        CCA_cum_count_allK = cumsum(CCA_component_count, 2);
        G_FWCCA_cum_count_allK = cumsum(G_FWCCA_component_count, 2);


        %% ====================================================
        % K-specific local-method analysis
        %
        % LG / GL are refitted because their CV-selected
        % bandwidth h depends on cumulative K.
        % ====================================================

        for id_k = 1:Kcomp_max

            fprintf('\n');
            fprintf('############################################################\n');
            fprintf('View = %s | p%d | cumulative K = %d\n', ...
                viewType, percentileValue, id_k);
            fprintf('############################################################\n');


            %% =================================================
            % K-specific CV-selected (h,r)
            %% =================================================

            LG_hr = LocalFWCCAoutput(id_p).LG_best_hAndr(id_k, :);
            GL_hr = LocalFWCCAoutput(id_p).GL_best_hAndr(id_k, :);

            LG_h = LG_hr(1);
            LG_r = LG_hr(2);
            GL_h = GL_hr(1);
            GL_r = GL_hr(2);

            fprintf('LG-FWCCA selected: h = %.3f, r = %.3f\n', LG_h, LG_r);
            fprintf('GL-FWCCA selected: h = %.3f, r = %.3f\n', GL_h, GL_r);


            %% =================================================
            % Construct K-specific local kernels
            %% =================================================

            K_LG = makeLocalKernel(T_train, LG_r, 'gaussian', LG_h, TR);
            K_GL = makeLocalKernel(T_train, GL_r, 'gaussian', GL_h, TR);


            %% =================================================
            % Refit LG-FWCCA
            %
            % LG = local kernel followed by global scaling
            %    = postscale
            %% =================================================

            res_LG = Run_FWCCA_Family( ...
                X1Train, X2Train, Wg_0, Wg_0, K_LG, K_LG, ...
                id_k, 100, 1e-5, 'postscale', ridgeTol, svdTol);


            %% =================================================
            % Refit GL-FWCCA
            %
            % GL = global scaling followed by local kernel
            %    = prescale
            %% =================================================

            res_GL = Run_FWCCA_Family( ...
                X1Train, X2Train, Wg_0, Wg_0, K_GL, K_GL, ...
                id_k, 100, 1e-5, 'prescale', ridgeTol, svdTol);


            %% =================================================
            % Project TEST data
            %
            % Output size: Nvox x id_k
            %% =================================================

            LG_FWCCAmaps_test = (XMask_test - res_LG.X1_localMean) * ...
                res_LG.Ak_raw(:, 1:id_k);

            GL_FWCCAmaps_test = (XMask_test - res_GL.X1_localMean) * ...
                res_GL.Ak_raw(:, 1:id_k);


            %% =================================================
            % Initialize cumulative local-method counts
            %% =================================================

            LG_FWCCA_cum_count = zeros(nSystems, 1);
            GL_FWCCA_cum_count = zeros(nSystems, 1);


            %% =================================================
            % Accumulate local-method counts over components 1:K
            %% =================================================

            for comp = 1:id_k

                fprintf('Accounting local component %d / %d\n', comp, id_k);


                %% ---------------------------------------------
                % LG-FWCCA
                %% ---------------------------------------------

                canonical_map = LG_FWCCAmaps_test(:, comp);

                [~, system_table] = CanonicalMap_to_FunctionalSystem( ...
                    canonical_map, idx_mask, vol_size, Vnii, atlas, ...
                    label_names, label_systems, systems_ordered, nTopVoxels, comp);

                LG_FWCCA_cum_count = LG_FWCCA_cum_count + system_table.Count;


                %% ---------------------------------------------
                % GL-FWCCA
                %% ---------------------------------------------

                canonical_map = GL_FWCCAmaps_test(:, comp);

                [~, system_table] = CanonicalMap_to_FunctionalSystem( ...
                    canonical_map, idx_mask, vol_size, Vnii, atlas, ...
                    label_names, label_systems, systems_ordered, nTopVoxels, comp);

                GL_FWCCA_cum_count = GL_FWCCA_cum_count + system_table.Count;

            end


            %% =================================================
            % Global cumulative counts for current K
            %% =================================================

            CCA_cum_count = CCA_cum_count_allK(:, id_k);
            G_FWCCA_cum_count = G_FWCCA_cum_count_allK(:, id_k);


            %% =================================================
            % Four-method comparison table
            %% =================================================

            comparison_table = table( ...
                string(systems_ordered(:)), ...
                CCA_cum_count, ...
                G_FWCCA_cum_count, ...
                LG_FWCCA_cum_count, ...
                GL_FWCCA_cum_count, ...
                'VariableNames', {'System', 'CCA', 'G_FWCCA', 'LG_FWCCA', 'GL_FWCCA'});


            %% =================================================
            % Display
            %% =================================================

            fprintf('\n');
            fprintf('Cumulative functional-system comparison\n');
            fprintf('Subject=%s | View=%s | K=%d | p%d\n', ...
                subj, viewType, id_k, percentileValue);

            disp(comparison_table);


            %% =================================================
            % Output workbook
            %
            % Results/
            %   MDD/
            %     sub-mdd02/
            %       K3/
            %         positive_MethodComparison.xlsx
            %         negative_MethodComparison.xlsx
            %% =================================================

            K_outputRoot = fullfile(outputRoot, sprintf('K%d', id_k));
            comparisonWorkbook = fullfile(K_outputRoot, ...
                sprintf('%s_MethodComparison.xlsx', viewType));


            %% =================================================
            % Current threshold -> one Excel sheet
            %
            % p90, p85, p80, p75, p70, p65
            %% =================================================

            comparisonSheet = sprintf('p%d', percentileValue);

            writetable(comparison_table, comparisonWorkbook, ...
                'Sheet', comparisonSheet);

            fprintf('Saved K=%d, %s, p%d\n', ...
                id_k, viewType, percentileValue);

        end

    end

end


%% ============================================================
% Finished
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Completed subject %s\n', subj);
fprintf('Results saved under:\n%s\n', outputRoot);
fprintf('============================================================\n');


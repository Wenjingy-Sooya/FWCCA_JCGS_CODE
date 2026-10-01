%% ============================================================
% Plot_ND_vs_MDD_KSpecific_FunctionalSystems.m
%
% Plot threshold-averaged functional-system counts for the Control (ND)
% and MDD groups using CCA, G-FWCCA, LG-FWCCA, and GL-FWCCA.
%
% Group comparisons are generated for:
%
%   K = 3, 4, 5
%
% separately for:
%
%   positive view
%   negative view
%
%
% INPUT
% -----
% Group-level threshold-averaged summaries:
%
%   Results/
%       ND/ND_Group_ThresholdAveraged_KSpecificSummary.xlsx
%       MDD/MDD_Group_ThresholdAveraged_KSpecificSummary.xlsx
%
% Each workbook contains the following sheets:
%
%   K3_positive
%   K3_negative
%   K4_positive
%   K4_negative
%   K5_positive
%   K5_negative
%
% Each sheet contains:
%
%   System
%   CCA_Mean
%   CCA_SD
%   G_FWCCA_Mean
%   G_FWCCA_SD
%   LG_FWCCA_Mean
%   LG_FWCCA_SD
%   GL_FWCCA_Mean
%   GL_FWCCA_SD
%
%
% FIGURE REPRESENTATION
% ---------------------
% Control (ND):
%   blue circle
%
% MDD:
%   pink/red triangle
%
% Each figure displays six functional systems:
%
%   Emotion Control
%   Visual
%   Auditory
%   Language/Cognition
%   Memory
%   Attention
%
% The x-axis compares:
%
%   CCA
%   G-FWCCA
%   LG-FWCCA
%   GL-FWCCA
%
%
% MANUSCRIPT FIGURES
% ------------------
% The K = 5 results are reported in the main manuscript:
%
%   Figure 11:
%   Average functional-system counts across thresholds for the control
%   and MDD groups for:
%       (a) positive view
%       (b) negative view
%   with K = 5.
%
% The K = 3 and K = 4 results are reported in the Supplementary Material:
%
%   Supplementary Figure 17:
%   Average functional-system counts across percentile thresholds for
%   the control and MDD groups with K = 3:
%       (a) positive view
%       (b) negative view.
%
%   Supplementary Figure 18:
%   Average functional-system counts across percentile thresholds for
%   the control and MDD groups with K = 4:
%       (a) positive view
%       (b) negative view.
%
%
% OUTPUT
% ------
% Figures are saved under:
%
%   Figures/
%
% Six figure files are generated:
%
%   K3 positive   -> Supplementary Figure 17(a)
%   K3 negative   -> Supplementary Figure 17(b)
%
%   K4 positive   -> Supplementary Figure 18(a)
%   K4 negative   -> Supplementary Figure 18(b)
%
%   K5 positive   -> Main-text Figure 11(a)
%   K5 negative   -> Main-text Figure 11(b)
%
% Each figure is saved in:
%
%   PDF (vector)
%   PNG (300 dpi)
%% ============================================================

clear;
clc;
close all;


%% ============================================================
% Evaluation paths
%% ============================================================

evaluationRoot = ...
    fileparts(mfilename('fullpath'));

resultsRoot = fullfile( ...
    evaluationRoot, ...
    'Results');

figureRoot = fullfile( ...
    evaluationRoot, ...
    'Figures');


if ~exist(resultsRoot, 'dir')

    error( ...
        'Results directory does not exist:\n%s', ...
        resultsRoot);

end


if ~exist(figureRoot, 'dir')

    mkdir(figureRoot);

end


%% ============================================================
% Input threshold-averaged group-summary files
%% ============================================================

ND_file = fullfile( ...
    resultsRoot, ...
    'ND', ...
    'ND_Group_ThresholdAveraged_KSpecificSummary.xlsx');

MDD_file = fullfile( ...
    resultsRoot, ...
    'MDD', ...
    'MDD_Group_ThresholdAveraged_KSpecificSummary.xlsx');


if ~exist(ND_file, 'file')

    error( ...
        'Cannot find ND summary file:\n%s', ...
        ND_file);

end


if ~exist(MDD_file, 'file')

    error( ...
        'Cannot find MDD summary file:\n%s', ...
        MDD_file);

end


fprintf('\n');
fprintf('============================================================\n');
fprintf('ND summary file:\n%s\n', ND_file);
fprintf('\n');
fprintf('MDD summary file:\n%s\n', MDD_file);
fprintf('============================================================\n');


%% ============================================================
% Numbers of cumulative canonical components
%% ============================================================

K_list = [3 4 5];

nK = numel(K_list);


%% ============================================================
% Views
%% ============================================================

view_list = { ...
    'positive', ...
    'negative'};

view_titles = { ...
    'Positive', ...
    'Negative'};

nViews = numel(view_list);


%% ============================================================
% Methods
%% ============================================================

methods = { ...
    'CCA', ...
    'G_FWCCA', ...
    'LG_FWCCA', ...
    'GL_FWCCA'};

method_labels = { ...
    'CCA', ...
    'G-FWCCA', ...
    'LG-FWCCA', ...
    'GL-FWCCA'};

nMethods = numel(methods);


%% ============================================================
% Functional systems
%
% Motor is excluded because its representation is negligible.
%% ============================================================

systems_plot = { ...
    'EmotionControl', ...
    'Visual', ...
    'Auditory', ...
    'LanguageCognition', ...
    'Memory', ...
    'Attention'};

system_titles = { ...
    'Emotion Control', ...
    'Visual', ...
    'Auditory', ...
    'Language/Cognition', ...
    'Memory', ...
    'Attention'};

nSystems = numel(systems_plot);


%% ============================================================
% Group colors and markers
%% ============================================================

ND_color = [ ...
    0.20 0.55 0.80];

MDD_color = [ ...
    0.90 0.35 0.40];

ND_marker = 'o';

MDD_marker = '^';


%% ============================================================
% Plot settings
%% ============================================================

fs_tick   = 22;
fs_title  = 22;
fs_ylabel = 22;
fs_main   = 22;
fs_legend = 22;

markerSize_ND  = 10;
markerSize_MDD = 10;

markerLineWidth = 1;

figure_width  = 14;
figure_height = 7.8;

% Slight horizontal offset so Control and MDD do not overlap
groupOffset = 0.15;


%% ============================================================
% Loop over K = 3, 4, 5
%% ============================================================

for kIdx = 1:nK

    K = K_list(kIdx);


    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('Processing K = %d\n', K);
    fprintf('============================================================\n');


    %% ========================================================
    % Loop over positive / negative views
    %% ========================================================

    for v = 1:nViews

        viewType  = view_list{v};
        viewTitle = view_titles{v};


        %% ----------------------------------------------------
        % Corresponding Excel sheet
        %
        % Example:
        %
        %   K3_positive
        %   K3_negative
        %% ----------------------------------------------------

        sheetName = sprintf( ...
            'K%d_%s', ...
            K, ...
            viewType);


        fprintf('\n');
        fprintf('------------------------------------------------------------\n');
        fprintf('K = %d | View = %s\n', ...
            K, ...
            viewTitle);
        fprintf('Reading sheet: %s\n', ...
            sheetName);
        fprintf('------------------------------------------------------------\n');


        %% ====================================================
        % Read group-summary tables
        %% ====================================================

        ND_tbl = readtable( ...
            ND_file, ...
            'Sheet', ...
            sheetName);

        MDD_tbl = readtable( ...
            MDD_file, ...
            'Sheet', ...
            sheetName);


        %% ====================================================
        % Create figure
        %% ====================================================

        fig = figure( ...
            'Name', ...
            sprintf( ...
                'ND_vs_MDD_K%d_%sView_FunctionalSystems', ...
                K, ...
                viewTitle), ...
            'Color', ...
            'w', ...
            'Units', ...
            'inches', ...
            'Position', ...
            [0.5 0.5 figure_width figure_height]);


        tl = tiledlayout( ...
            2, ...
            3, ...
            'TileSpacing', ...
            'loose', ...
            'Padding', ...
            'compact');


        %% ====================================================
        % One subplot per functional system
        %% ====================================================

        for s = 1:nSystems

            ax = nexttile;

            hold(ax, 'on');


            systemName = ...
                systems_plot{s};


            %% ------------------------------------------------
            % Locate system row
            %% ------------------------------------------------

            ND_row = strcmpi( ...
                strtrim(string(ND_tbl.System)), ...
                systemName);

            MDD_row = strcmpi( ...
                strtrim(string(MDD_tbl.System)), ...
                systemName);


            if ~any(ND_row)

                warning( ...
                    'System %s not found in ND table.', ...
                    systemName);

                continue;

            end


            if ~any(MDD_row)

                warning( ...
                    'System %s not found in MDD table.', ...
                    systemName);

                continue;

            end


            %% ------------------------------------------------
            % Extract method-specific means
            %% ------------------------------------------------

            ND_mean = ...
                zeros( ...
                    nMethods, ...
                    1);

            MDD_mean = ...
                zeros( ...
                    nMethods, ...
                    1);


            for m = 1:nMethods

                meanVar = sprintf( ...
                    '%s_Mean', ...
                    methods{m});


                if ~ismember( ...
                        meanVar, ...
                        ND_tbl.Properties.VariableNames)

                    error( ...
                        'Variable %s not found in ND sheet %s.', ...
                        meanVar, ...
                        sheetName);

                end


                if ~ismember( ...
                        meanVar, ...
                        MDD_tbl.Properties.VariableNames)

                    error( ...
                        'Variable %s not found in MDD sheet %s.', ...
                        meanVar, ...
                        sheetName);

                end


                ND_mean(m) = ...
                    ND_tbl{ ...
                        ND_row, ...
                        meanVar};

                MDD_mean(m) = ...
                    MDD_tbl{ ...
                        MDD_row, ...
                        meanVar};

            end


            %% =================================================
            % x locations
            %% =================================================

            x = 1:nMethods;

            xControl = ...
                x - groupOffset;

            xMDD = ...
                x + groupOffset;


            %% =================================================
            % Control markers
            %% =================================================

            plot( ...
                ax, ...
                xControl, ...
                ND_mean, ...
                ND_marker, ...
                'LineStyle', ...
                'none', ...
                'MarkerSize', ...
                markerSize_ND, ...
                'MarkerFaceColor', ...
                ND_color, ...
                'MarkerEdgeColor', ...
                'black', ...
                'LineWidth', ...
                markerLineWidth);


            %% =================================================
            % MDD markers
            %% =================================================

            plot( ...
                ax, ...
                xMDD, ...
                MDD_mean, ...
                MDD_marker, ...
                'LineStyle', ...
                'none', ...
                'MarkerSize', ...
                markerSize_MDD, ...
                'MarkerFaceColor', ...
                MDD_color, ...
                'MarkerEdgeColor', ...
                'black', ...
                'LineWidth', ...
                markerLineWidth);


            %% =================================================
            % Axis formatting
            %% =================================================

            set( ...
                ax, ...
                'XTick', ...
                1:nMethods, ...
                'XTickLabel', ...
                method_labels, ...
                'FontSize', ...
                fs_tick, ...
                'LineWidth', ...
                1.1);


            xtickangle( ...
                ax, ...
                25);


            %% ------------------------------------------------
            % Subplot title
            %% ------------------------------------------------

            title( ...
                ax, ...
                system_titles{s}, ...
                'FontSize', ...
                fs_title, ...
                'FontWeight', ...
                'bold');


            %% ------------------------------------------------
            % y-axis
            %% ------------------------------------------------

            ylabel( ...
                ax, ...
                'Average count', ...
                'FontSize', ...
                fs_ylabel);


            %% ------------------------------------------------
            % x-limits
            %% ------------------------------------------------

            xlim( ...
                ax, ...
                [0.5, nMethods + 0.5]);


            %% ------------------------------------------------
            % Use the same y-axis range for all systems
            %% ------------------------------------------------

            ylim( ...
                ax, ...
                [0 16]);

            yticks( ...
                ax, ...
                0:4:16);


            %% ------------------------------------------------
            % Appearance
            %% ------------------------------------------------

            box(ax, 'on');

            grid(ax, 'off');

            hold(ax, 'off');

        end


        %{
        %% ====================================================
        % Main title
        %
        % Kept disabled to match the original figure style.
        %% ====================================================

        title( ...
            tl, ...
            sprintf( ...
                'Control vs MDD (%s View)', ...
                viewTitle), ...
            'FontSize', ...
            fs_main, ...
            'FontWeight', ...
            'bold');

        %}


        %% ====================================================
        % Shared legend
        %
        % Use dummy marker handles so the legend is independent
        % of the plotted data in the six panels.
        %% ====================================================

        hold(ax, 'on');


        hLegendControl = plot( ...
            ax, ...
            NaN, ...
            NaN, ...
            ND_marker, ...
            'LineStyle', ...
            'none', ...
            'MarkerSize', ...
            markerSize_ND, ...
            'MarkerFaceColor', ...
            ND_color, ...
            'MarkerEdgeColor', ...
            'black', ...
            'LineWidth', ...
            markerLineWidth);


        hLegendMDD = plot( ...
            ax, ...
            NaN, ...
            NaN, ...
            MDD_marker, ...
            'LineStyle', ...
            'none', ...
            'MarkerSize', ...
            markerSize_MDD, ...
            'MarkerFaceColor', ...
            MDD_color, ...
            'MarkerEdgeColor', ...
            'black', ...
            'LineWidth', ...
            markerLineWidth);


        lgd = legend( ...
            [hLegendControl, hLegendMDD], ...
            { ...
                'Control', ...
                'MDD'}, ...
            'Orientation', ...
            'horizontal', ...
            'FontSize', ...
            fs_legend, ...
            'Box', ...
            'on');


        lgd.Layout.Tile = ...
            'south';


        hold(ax, 'off');


        %% ====================================================
        % Save figure
        %% ====================================================

        safeName = sprintf( ...
            ['Control_vs_MDD_K%d_%sView_' ...
             'ThresholdAveraged_FunctionalSystems_Markers'], ...
            K, ...
            viewTitle);


        %% ----------------------------------------------------
        % PDF
        %% ----------------------------------------------------

        pdfFile = fullfile( ...
            figureRoot, ...
            [safeName '.pdf']);


        exportgraphics( ...
            fig, ...
            pdfFile, ...
            'ContentType', ...
            'vector');


        %% ----------------------------------------------------
        % PNG
        %% ----------------------------------------------------

        pngFile = fullfile( ...
            figureRoot, ...
            [safeName '.png']);


        exportgraphics( ...
            fig, ...
            pngFile, ...
            'Resolution', ...
            300);


        fprintf('\n');
        fprintf('============================================================\n');
        fprintf( ...
            'K = %d | %s view completed.\n', ...
            K, ...
            viewTitle);

        fprintf( ...
            'PDF saved to:\n%s\n', ...
            pdfFile);

        fprintf( ...
            'PNG saved to:\n%s\n', ...
            pngFile);

        fprintf('============================================================\n');


    end

end


%% ============================================================
% Finished
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('Control vs MDD marker comparison completed.\n');
fprintf('K values: ');
fprintf('%d ', K_list);
fprintf('\n');
fprintf('Figures saved under:\n%s\n', ...
    figureRoot);
fprintf('============================================================\n');
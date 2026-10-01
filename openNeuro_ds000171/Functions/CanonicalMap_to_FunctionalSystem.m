function [region_table, system_table] = ...
    CanonicalMap_to_FunctionalSystem( ...
        canonical_map, ...
        idx_mask, ...
        vol_size, ...
        Vnii, ...
        atlas, ...
        label_names, ...
        label_systems, ...
        systems_ordered,...
        nTopVoxels, ...
        compIndex)

% ============================================================
% CanonicalMap_to_FunctionalSystem
%
% Map one spatial canonical component to anatomical regions
% and functional systems.
%
% INPUTS
%
% canonical_map : [Nvox x 1]
%                 Spatial canonical map within the analysis mask.
%
% idx_mask      : Linear indices of analysis-mask voxels in the
%                 full 3-D image volume.
%
% vol_size      : Size of the full 3-D analysis volume.
%
% Vnii          : SPM NIfTI header used for voxel-to-MNI
%                 coordinate transformation.
%
% atlas         : SPM atlas object.
%
% label_names   : Anatomical-region names used in the
%                 region-to-functional-system mapping.
%
% label_systems : Corresponding functional-system labels.
%                 Multiple systems may be separated by ';'.
%
% topProportion : Proportion of voxels retained from the absolute
%                 canonical map, e.g. 0.05 for top 5%.
%
% compIndex     : Canonical-component index used for output tables.
%
%
% OUTPUTS
%
% region_table  : Table with columns
%                 Comp, Region, Count, Systems
%
% system_table  : Table with columns
%                 Comp, System, Count
% ============================================================


%% ============================================================
% Basic checks
%% ============================================================

canonical_map = canonical_map(:);

if numel(canonical_map) ~= numel(idx_mask)

    error(['canonical_map must have the same number of elements ' ...
           'as idx_mask.']);

end


if nargin < 9 || isempty(compIndex)

    compIndex = 1;

end



%% ============================================================
% Absolute canonical map
%% ============================================================

comp = abs(canonical_map);

max_comp = max(comp);

if max_comp > eps

    comp_mask = comp / max_comp;

else

    comp_mask = comp;

end


%% ============================================================
% Rank voxels and retain the top proportion
%% ============================================================

[~, sort_idx] = sort(comp_mask, 'descend');

nTopVoxels = min(nTopVoxels, numel(sort_idx));
selected_idx = sort_idx(1:nTopVoxels);


%% ============================================================
% Query SPM atlas
%% ============================================================

all_region_names = {};

for k = 1:nTopVoxels

    %% --------------------------------------------------------
    % Mask-space index -> full-volume linear index
    %% --------------------------------------------------------

    voxel_idx = ...
        idx_mask(selected_idx(k));

    %% --------------------------------------------------------
    % Full-volume linear index -> voxel coordinates
    %% --------------------------------------------------------

    [x,y,z] = ...
        ind2sub( ...
            vol_size, ...
            voxel_idx);

    %% --------------------------------------------------------
    % Voxel coordinates -> MNI coordinates
    %% --------------------------------------------------------

    mni_coords = ...
        Vnii(1).mat * ...
        [x; y; z; 1];

    %% --------------------------------------------------------
    % Query Neuromorphometrics atlas
    %% --------------------------------------------------------

    xY = struct( ...
        'xyz', ...
        mni_coords(1:3), ...
        'def', ...
        'sphere', ...
        'spec', ...
        1);

    region = ...
        spm_atlas( ...
            'query', ...
            atlas, ...
            xY);

    if isempty(region)

        continue;

    end

    all_region_names{end+1} = ...
        strjoin( ...
            region, ...
            ', ');

end


%% ============================================================
% Clean atlas region labels
%% ============================================================

if isempty(all_region_names)

    flat_names = {};

else

    flat_names = ...
        strsplit( ...
            strjoin( ...
                all_region_names, ...
                ', '), ...
            ', ');

    flat_names = ...
        strtrim(flat_names);

    flat_names = ...
        flat_names( ...
            ~cellfun( ...
                @isempty, ...
                flat_names));

    flat_names = ...
        flat_names( ...
            ~contains( ...
                flat_names, ...
                'Unknown'));

end


%% ============================================================
% Count anatomical-region frequencies
%% ============================================================

if isempty(flat_names)

    regions_sorted = {};
    region_counts_sorted = [];

else

    [uniq_regions,~,ic] = ...
        unique(flat_names);

    region_counts = ...
        accumarray( ...
            ic(:), ...
            1);

    [region_counts_sorted,idx_sort] = ...
        sort( ...
            region_counts, ...
            'descend');

    regions_sorted = ...
        uniq_regions(idx_sort);

end


%% ============================================================
% Map anatomical regions to functional systems
%% ============================================================

system_freq_map = ...
    containers.Map( ...
        'KeyType', ...
        'char', ...
        'ValueType', ...
        'double');

region_rows = {};

for i = 1:numel(regions_sorted)

    region_name = ...
        regions_sorted{i};

    region_count = ...
        region_counts_sorted(i);

    match_idx = ...
        strcmpi( ...
            label_names, ...
            region_name);

    matched_systems = {};

    %% --------------------------------------------------------
    % Region found in functional-system mapping
    %% --------------------------------------------------------

    if any(match_idx)

        matched_row = ...
            find( ...
                match_idx, ...
                1);

        systems_raw = ...
            label_systems{matched_row};

        matched_systems = ...
            strtrim( ...
                strsplit( ...
                    systems_raw, ...
                    ';'));

        for s = 1:numel(matched_systems)

            sys = ...
                matched_systems{s};

            if isempty(sys)

                continue;

            end

            if isKey( ...
                    system_freq_map, ...
                    sys)

                system_freq_map(sys) = ...
                    system_freq_map(sys) ...
                    + region_count;

            else

                system_freq_map(sys) = ...
                    region_count;

            end

        end

    else

        matched_systems = ...
            {'Unknown'};

    end

    %% --------------------------------------------------------
    % Store region-level result
    %% --------------------------------------------------------

    region_rows(end+1,:) = { ...
        compIndex, ...
        region_name, ...
        region_count, ...
        strjoin( ...
            matched_systems, ...
            ', ')};

end


%% ============================================================
% Region-frequency table
%% ============================================================

if isempty(region_rows)

    region_rows = { ...
        compIndex, ...
        'None', ...
        0, ...
        'Unknown'};

end

region_table = ...
    cell2table( ...
        region_rows, ...
        'VariableNames', { ...
            'Comp', ...
            'Region', ...
            'Count', ...
            'Systems'});


%% ============================================================
% Functional-system frequency table
%
% The output order is fixed according to systems_ordered.
% Systems not detected in the current canonical component
% are assigned Count = 0.
%% ============================================================

nSystems = numel(systems_ordered);

system_counts = zeros(nSystems,1);

for s = 1:nSystems

    sys = systems_ordered{s};

    if isKey(system_freq_map, sys)

        system_counts(s) = ...
            system_freq_map(sys);

    else

        system_counts(s) = 0;

    end

end

system_table = table( ...
    repmat(compIndex,nSystems,1), ...
    string(systems_ordered(:)), ...
    system_counts, ...
    'VariableNames', { ...
        'Comp', ...
        'System', ...
        'Count'});

end
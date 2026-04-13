%% fig4_phylogenetic_tree.m
% Core genome phylogenetic tree colored by protection score
% Rerooted at VPI then VPI pruned — shows ST1 strain variation only
% Run with MATLAB working directory set to data/raw/genomics/

%% Load protection scores
tbl         = readtable('binary_gene_pre_abs.csv');
strainNames = string(tbl{:, 1});
strainNames = strrep(strainNames, 'ST1-', 'ST1.');
protection  = tbl{:, end};

% Color map: pure green (low protection) to pure blue (high protection)
protCmap = [zeros(256,1), linspace(0.8,0,256)', linspace(0,1,256)'];
minP = min(protection);
maxP = max(protection);

%% Load strain name mapping
names = readtable('strain_names_mapping.csv', 'ReadVariableNames', false);

%% Read newick and rename strains
fid = fopen('VPI Outgroup_tree.nwk', 'r');
newick_str = fread(fid, '*char')';
fclose(fid);

% Replace numeric IDs with ST1.XX names
for i = 1:height(names)
    for decimals = [3 4 5 6]
        old_name = sprintf(['%.' num2str(decimals) 'f'], names.Var1(i));
        new_name = strrep(char(string(names.Var2(i))), 'ST1-', 'ST1.');
        newick_str = strrep(newick_str, old_name, new_name);
    end
end

% Fix VPI naming
newick_str = strrep(newick_str, '1496.34290', 'VPI');
newick_str = strrep(newick_str, '1496.343',   'VPI');
newick_str = strrep(newick_str, 'VPI0',        'VPI');

% Write renamed tree
fid = fopen('tree_renamed.nwk', 'w');
fprintf(fid, '%s', newick_str);
fclose(fid);

%% Read tree
tree = phytreeread('tree_renamed.nwk');
leaf_names = get(tree, 'LeafNames');

fprintf('\nLeaf names after renaming:\n');
disp(leaf_names);

%% Reroot at VPI
vpi_idx = find(strcmp(leaf_names, 'VPI'));
if isempty(vpi_idx)
    fprintf('WARNING: VPI not found — check renaming\n');
else
    tree = reroot(tree, vpi_idx);
    fprintf('Rerooted at VPI successfully\n');

    % Prune VPI so only ST1 strains appear in figure
    leaf_names = get(tree, 'LeafNames');
    vpi_idx    = find(strcmp(leaf_names, 'VPI'));
    tree       = prune(tree, vpi_idx);
    fprintf('VPI pruned — plotting ST1 strains only\n');
end

%% Plot tree
figure('Position', [100 100 900 700]);
h = plot(tree, 'Type', 'square');
set(findall(gca, 'Type', 'text'), 'FontSize', 12, 'FontName', 'Arial');
set(findall(gca, 'Type', 'line'), 'LineWidth', 2);
set(gca, 'FontSize', 15);
set(h.BranchDots, 'MarkerSize', 8);
set(h.LeafDots,   'MarkerSize', 8);

% Give more room for strain labels on right
ax = gca;
ax.Position(3) = ax.Position(3) * 0.75;

%% Color leaf nodes by protection score
text_handles = findall(gca, 'Type', 'text');

hold on;
matched = 0;
for i = 1:length(text_handles)
    label = get(text_handles(i), 'String');

    % Match label to strain — both in ST1.XX format
    idx = find(strcmp(strainNames, label));

    if isempty(idx)
        idx = find(strcmp(strrep(strainNames, 'ST1.', 'ST1-'), label));
    end

    if ~isempty(idx)
        matched = matched + 1;
        prot_val = protection(idx(1));
        colorIdx = round((prot_val - minP) / (maxP - minP) * 255) + 1;
        colorIdx = max(1, min(256, colorIdx));

        set(text_handles(i), 'Color', 'k');
        pos = get(text_handles(i), 'Position');
        scatter(pos(1), pos(2), 120, protCmap(colorIdx,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end
end
hold off;

fprintf('Colored %d out of %d leaf labels\n', matched, length(text_handles));

%% Colorbar and formatting
colormap(protCmap);
cb = colorbar;
cb.Label.String   = 'Protection score';
cb.Label.FontSize = 14;
cb.FontSize       = 12;
clim([minP maxP]);
title('Core genome phylogeny colored by protection score', 'FontSize', 14);
set(gca, 'FontSize', 15, 'FontName', 'Arial');
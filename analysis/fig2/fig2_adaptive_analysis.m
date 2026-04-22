

%% Load all cell data
mice = [1 6 7 11 12 13];

allMice = [];
for i = 1:length(mice)
    allCellsMX = readtable(sprintf('ks10_adaptive_csv_files/all_adaptive/all_adaptive_Specimen_001_%d.csv', mice(i)));
    allCellsMX.mouse(:) = mice(i);
    if i == 1
        allMice = allCellsMX;
    else
        allMice = [allMice; allCellsMX];
    end
end

%% Add gate flags
allMiceWithGates = allMice;
cellTypes = {'B_cells' 'CD19negTCRbneg' 'NK_cells' 'T_cells' 'CD4pos' 'CD8pos'};

for j = 1:length(cellTypes)
    allCellsOfType = [];
    dirName = sprintf('ks10_adaptive_csv_files/%s/*.csv', cellTypes{j});
    d = dir(dirName);
    for i = 1:length(d)
        allCellsMX = readtable(sprintf('ks10_adaptive_csv_files/%s/%s', cellTypes{j}, d(i).name));
        aux = strsplit(d(i).name, '_001_');
        aux = strsplit(aux{2}, '_');
        mouse = str2double(aux{1});
        allCellsMX.mouse(:) = mouse;
        if i == 1
            allCellsOfType = allCellsMX;
        else
            allCellsOfType = [allCellsOfType; allCellsMX];
        end
    end
    binaryGate = ismember(allMice, allCellsOfType);
    allMiceWithGates = addvars(allMiceWithGates, binaryGate);
    allMiceWithGates.Properties.VariableNames{end} = cellTypes{j};
end

%% Remove events that did not pass any gate
didNotPassAnyGate = sum(allMiceWithGates{:, 17:22}, 2) == 0;
allGatedCell = allMiceWithGates(didNotPassAnyGate == 0, :);

%% Find cells within fluorescence range of gated cells
fluor = allGatedCell{:, 7:14};
allFlowMetrics = allGatedCell{:, 1:14};

maxValuesAmongGated = max(allFlowMetrics);
minValuesAmongGated = min(allFlowMetrics);
allFlowMetricsAllCells = allMiceWithGates{:, 1:14};
cellsWithinRanges = all(allFlowMetricsAllCells >= minValuesAmongGated & ...
    allFlowMetricsAllCells <= maxValuesAmongGated, 2);
allCellsWithinRange = allMiceWithGates(cellsWithinRanges, :);

fluor = allCellsWithinRange{:, 7:14};
allFlowMetrics = allCellsWithinRange{:, 1:14};

%% Define mouse groups
UI_mice = [1, 6, 7];
Avirulent_mice = [11, 12, 13];

mouseGroup = repmat({'Other'}, height(allCellsWithinRange), 1);
mouseGroup(ismember(allCellsWithinRange.mouse, UI_mice)) = {'UI'};
mouseGroup(ismember(allCellsWithinRange.mouse, Avirulent_mice)) = {'Avirulent'};
mouseGroup = categorical(mouseGroup);
prettyTypeNames = {'B cells', 'CD19^{-}TCR\beta^{-}', 'NK cells', ...
    'T cells', 'CD4^{+}', 'CD8^{+}'};

%% Generate categorical cell type labels
cellTypeVars = allCellsWithinRange.Properties.VariableNames(17:22);

cellLabels = repmat({'Other'}, height(allCellsWithinRange), 1);
for i = 1:length(cellTypeVars)
    idx = allCellsWithinRange{:, cellTypeVars{i}} == 1;
    cellLabels(idx) = cellTypeVars(i);
end
cellLabels = categorical(cellLabels);

%% PCA
[coeff, tscore_pca, latent, ~, explained, mu] = pca(fluor);

figure;
hold on;

namedCats_pca = categories(cellLabels);
namedCats_pca(strcmp(namedCats_pca, 'Other')) = [];
colors_pca = lines(numel(namedCats_pca));

scatter(tscore_pca(cellLabels=='Other', 1), tscore_pca(cellLabels=='Other', 2), ...
    10, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.3);

for i = 1:length(namedCats_pca)
    idx = cellLabels == categorical(namedCats_pca(i));
    scatter(tscore_pca(idx,1), tscore_pca(idx,2), 15, colors_pca(i,:), ...
        'filled', 'MarkerFaceAlpha', 0.9);
end

xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('PCA: Adaptive immune cell types');

legendLabels_pca = [{'Other'}; namedCats_pca];
for i = 1:length(cellTypeVars)
    legendLabels_pca(strcmp(legendLabels_pca, cellTypeVars{i})) = prettyTypeNames(i);
end
lgd = legend(legendLabels_pca, 'Location', 'bestoutside');
set(lgd, 'Interpreter', 'tex');
set(lgd, 'FontSize', 14);
set(lgd, 'ItemTokenSize', [50, 18]);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
hold off;

%% PCA colored by infection group
figure;
hold on;
scatter(tscore_pca(mouseGroup=='UI', 1), tscore_pca(mouseGroup=='UI', 2), ...
    30, [0.2 0.6 1], 'filled', 'MarkerFaceAlpha', 0.7);
scatter(tscore_pca(mouseGroup=='Avirulent', 1), tscore_pca(mouseGroup=='Avirulent', 2), ...
    30, [0.6 0.8 0.3], 'filled', 'MarkerFaceAlpha', 0.7);
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('PCA: Infection group');
legend({'Uninfected (UI)', 'ST1.75 colonized'}, 'Location', 'bestoutside');
set(gca, 'FontName', 'Arial', 'FontSize', 20);
hold off;

%% tSNE — fixed seed for reproducibility
rng(42);
[tscore] = tsne(allFlowMetrics, 'Standardize', 1);

if length(cellLabels) ~= size(tscore, 1)
    error('Mismatch: %d labels vs %d tSNE points', length(cellLabels), size(tscore, 1));
end

%% tSNE colored by cell type
figure;
hold on;

otherIdx = cellLabels == categorical({'Other'});
scatter(tscore(otherIdx,1), tscore(otherIdx,2), 10, ...
    [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.3);

colors = lines(numel(categories(cellLabels)) - 1);
namedCats = categories(cellLabels);
namedCats(strcmp(namedCats, 'Other')) = [];

for i = 1:length(namedCats)
    idx = cellLabels == categorical(namedCats(i));
    scatter(tscore(idx,1), tscore(idx,2), 15, colors(i,:), ...
        'filled', 'MarkerFaceAlpha', 0.9);
end

xlabel('tSNE 1'); ylabel('tSNE 2');
title('tSNE: Adaptive immune cell types');

legendLabels = [{'Other'}; namedCats];
for i = 1:length(cellTypeVars)
    legendLabels(strcmp(legendLabels, cellTypeVars{i})) = prettyTypeNames(i);
end

lgd = legend(legendLabels, 'Location', 'bestoutside');
set(lgd, 'Interpreter', 'tex');
set(lgd, 'FontSize', 14);
set(lgd, 'ItemTokenSize', [50, 18]);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
hold off;

%% tSNE colored by infection group
figure;
hold on;
scatter(tscore(mouseGroup=='UI', 1), tscore(mouseGroup=='UI', 2), ...
    15, [0.2 0.6 1], 'filled', 'MarkerFaceAlpha', 0.7);
scatter(tscore(mouseGroup=='Avirulent', 1), tscore(mouseGroup=='Avirulent', 2), ...
    15, [0.6 0.8 0.3], 'filled', 'MarkerFaceAlpha', 0.7);
xlabel('tSNE 1'); ylabel('tSNE 2');
title('tSNE: Infection group');
legend({'Uninfected (UI)', 'ST1.75 colonized'}, 'Location', 'bestoutside');
set(gca, 'FontName', 'Arial', 'FontSize', 20);
hold off;

%% Compute cell type fractions per mouse
cellTypeNames = allCellsWithinRange.Properties.VariableNames(17:22);
mouseIDs = unique(allCellsWithinRange.mouse);
countCells = table(mouseIDs, 'VariableNames', {'mouse'});

for i = 1:length(mouseIDs)
    mouse = mouseIDs(i);
    dataMouse = allCellsWithinRange(allCellsWithinRange.mouse == mouse, :);
    totalEvents = height(dataMouse);
    fractions = zeros(1, length(cellTypeNames));
    for j = 1:length(cellTypeNames)
        fractions(j) = sum(dataMouse{:, cellTypeNames{j}}) / totalEvents;
    end
    countCells{i, 2:(length(cellTypeNames)+1)} = fractions;
end
countCells.Properties.VariableNames(2:end) = cellTypeNames;

%% Wilcoxon rank-sum tests
UI = [1 6 7];
Avirulent = [11 12 13];

data = countCells{:, 2:end};
mouseIDs = countCells.mouse;

UI_data = data(ismember(mouseIDs, UI), :);
Avirulent_data = data(ismember(mouseIDs, Avirulent), :);

UI_mean = mean(UI_data, 1);
Avirulent_mean = mean(Avirulent_data, 1);
UI_std = std(UI_data, 0, 1);
Avirulent_std = std(Avirulent_data, 0, 1);

groupMeans = [UI_mean; Avirulent_mean];
nCellTypes = length(cellTypeNames);

pvals = nan(1, nCellTypes);
for j = 1:nCellTypes
    pvals(j) = ranksum(UI_data(:,j), Avirulent_data(:,j));
end

fprintf('\nWilcoxon rank-sum p-values (n=3 per group):\n');
for j = 1:nCellTypes
    fprintf('  %s: p = %.4f\n', cellTypeNames{j}, pvals(j));
end

%% Bar plot 

figure;
hold on;

groupColors = [0.4 0.6 0.9; 0.9 0.5 0.3];
nGroups = 2;
groupWidth = 0.7;
barWidth = groupWidth / nGroups;

for i = 1:nGroups
    xPos = (1:nCellTypes) - groupWidth/2 + (i-0.5)*barWidth;
    bar(xPos, groupMeans(i,:), barWidth, 'FaceColor', groupColors(i,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.85);
end

for i = 1:nGroups
    xPos = (1:nCellTypes) - groupWidth/2 + (i-0.5)*barWidth;
    if i == 1
        groupData = UI_data;
        groupStd = UI_std;
    else
        groupData = Avirulent_data;
        groupStd = Avirulent_std;
    end
    errorbar(xPos, groupMeans(i,:), groupStd, 'k.', 'LineWidth', 1.5, 'CapSize', 5);
    for j = 1:nCellTypes
        scatter(repmat(xPos(j), size(groupData,1), 1), groupData(:,j), ...
            40, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    end
end

for j = 1:nCellTypes
    x1 = j - groupWidth/2 + 0.5*barWidth;
    x2 = j - groupWidth/2 + 1.5*barWidth;
    xmid = (x1 + x2) / 2;
    ymax = max(groupMeans(:,j)) + max(UI_std(j), Avirulent_std(j)) + 0.005;
    if pvals(j) < 0.001
        plot([x1 x2], [ymax ymax], 'k-', 'LineWidth', 1);
        text(xmid, ymax + 0.002, '***', 'FontSize', 16, ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    elseif pvals(j) < 0.01
        plot([x1 x2], [ymax ymax], 'k-', 'LineWidth', 1);
        text(xmid, ymax + 0.002, '**', 'FontSize', 16, ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    elseif pvals(j) < 0.05
        plot([x1 x2], [ymax ymax], 'k-', 'LineWidth', 1);
        text(xmid, ymax + 0.002, '*', 'FontSize', 16, ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

set(gca, 'XTick', 1:nCellTypes);
set(gca, 'XTickLabel', prettyTypeNames);
set(gca, 'XTickLabelRotation', 45);
set(gca, 'TickLabelInterpreter', 'tex');
legend({'Uninfected (UI)', 'ST1.75 colonized'}, 'Location', 'northeast');
set(gca, 'FontName', 'Arial', 'FontSize', 20);
xlabel('Adaptive immune cell type');
ylabel('Fraction of gated cells');
box on; grid on;
hold off;

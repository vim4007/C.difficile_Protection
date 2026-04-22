
%% Load data
tbl = readtable('/Users/vmishra/C.difficile_Protection/data/Genomics/binary_gene_pre_abs.csv');
strainNames = string(tbl{:, 1});
varNames    = tbl.Properties.VariableNames;
geneNames   = varNames(2:end-2);
data        = double(tbl{:, 2:end-2});
virulence   = tbl{:, end-1};
protection  = tbl{:, end};

nStrains = height(tbl);
nGenes   = length(geneNames);

[protection_sorted, sortIdx] = sort(protection, 'ascend');
strains_sorted   = strainNames(sortIdx);
virulence_sorted = virulence(sortIdx);
data_sorted      = data(sortIdx, :);
coreMask      = all(data == 1, 1);
accessoryMask = ~coreMask;
nCore         = sum(coreMask);
nAccessory    = sum(accessoryMask);

protCmap = [zeros(256,1), linspace(0.8,0,256)', linspace(0,1,256)'];
minP = min(protection);
maxP = max(protection);

%% Figure 1 — Pan-genome pie chart
figure('Position', [100 100 400 400]);
pie([nCore nAccessory]);
colormap([0 0.45 0.74; 0.9 0.5 0.3]);
legend({'Core genome', 'Accessory genome'}, 'Location', 'southoutside', 'FontSize', 20);
title(sprintf('Pan-genome: %d CDS', nGenes), 'FontSize', 20);

%% Figure 2 — PCA of accessory gene matrix colored by protection score
accessoryData          = data(:, accessoryMask);
accessoryData_centered = accessoryData - mean(accessoryData, 1);
[~, score, ~, ~, explained] = pca(accessoryData_centered);

figure('Position', [100 100 600 500]);
hold on;

for i = 1:nStrains
    colorIdx = round((protection(i) - minP) / (maxP - minP) * 255) + 1;
    scatter(score(i,1), score(i,2), 250, protCmap(colorIdx,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

for i = 1:nStrains
    text(score(i,1)+0.3, score(i,2)+0.3, ...
    strrep(strainNames{i}, 'ST1-', 'ST1.'), ...
    'FontSize', 25, 'Interpreter', 'none');
end

colormap(protCmap);
cb = colorbar;
cb.Label.String = 'Protection score';
clim([minP maxP]);
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));

set(gca, 'FontName', 'Arial', 'FontSize', 25);
box on; hold off;

%% Compute gene correlations with protection score

correlations = zeros(nGenes, 1);
pvals_gene   = ones(nGenes, 1);

for j = 1:nGenes
    if any(data(:,j) ~= data(1,j))
        [r, p] = corr(data(:,j), protection, 'Type', 'Spearman');
        correlations(j) = r;
        pvals_gene(j)   = p;
    end
end

%% Figure 4 — Waterfall plot top/bottom 10 gene correlations
bonferroni_threshold = 0.05 / nGenes;
t_thresh = tinv(1 - bonferroni_threshold/2, nStrains-2);
r_thresh = t_thresh / sqrt(t_thresh^2 + nStrains - 2);

validMask         = correlations ~= 0;
validCorr         = correlations(validMask);
validPvals        = pvals_gene(validMask);
validNames        = geneNames(validMask);

[validCorr_sorted, wSortIdx] = sort(validCorr, 'descend');
validPvals_sorted = validPvals(wSortIdx);
validNames_sorted = validNames(wSortIdx);
nValid            = length(validCorr_sorted);
plot_idx   = [1:10, (nValid-9):nValid];
plot_corr  = validCorr_sorted(plot_idx);
plot_pvals = validPvals_sorted(plot_idx);
plot_names = validNames_sorted(plot_idx);

% Colors
barColors = zeros(length(plot_idx), 3);
for i = 1:length(plot_idx)
    if plot_corr(i) >= 0
        barColors(i,:) = [0 0.45 0.74];
    else
        barColors(i,:) = [0.47 0.67 0.19];
    end
end

figure('Position', [100 100 950 520]);
hold on;

for i = 1:length(plot_idx)
    bar(i, plot_corr(i), 'FaceColor', barColors(i,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.85);
end

yline(r_thresh,  '--k', 'LineWidth', 1.5, ...
    'Label', sprintf('Bonferroni p < %.1e', bonferroni_threshold), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
yline(-r_thresh, '--k', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.5);

xline(10.5, 'k:', 'LineWidth', 1.5);

cleanNames = strrep(plot_names, '_', '\_');
set(gca, 'XTick', 1:length(plot_idx));
set(gca, 'XTickLabel', cleanNames);
set(gca, 'XTickLabelRotation', 45);
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'FontName', 'Arial', 'FontSize', 14);
xlabel('Gene (top 10 positive | bottom 10 negative correlation with protection)');
ylabel('Spearman r with protection score');
xlim([0 length(plot_idx)+1]);
box on; grid on; hold off;


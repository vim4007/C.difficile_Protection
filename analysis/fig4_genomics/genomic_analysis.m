%% fig4_genomics_analysis.m
% Genomics analysis — pan-genome, PCA, clcA1, waterfall plot
% Run with MATLAB working directory set to data/raw/genomics/

%% Load data
tbl = readtable('binary_gene_pre_abs.csv');
strainNames = string(tbl{:, 1});
varNames    = tbl.Properties.VariableNames;
geneNames   = varNames(2:end-2);
data        = double(tbl{:, 2:end-2});
virulence   = tbl{:, end-1};
protection  = tbl{:, end};

nStrains = height(tbl);
nGenes   = length(geneNames);

% Sort strains by protection score
[protection_sorted, sortIdx] = sort(protection, 'ascend');
strains_sorted   = strainNames(sortIdx);
virulence_sorted = virulence(sortIdx);
data_sorted      = data(sortIdx, :);

% Core and accessory
coreMask      = all(data == 1, 1);
accessoryMask = ~coreMask;
nCore         = sum(coreMask);
nAccessory    = sum(accessoryMask);

fprintf('\n=== Pan-genome summary ===\n');
fprintf('Total strains:     %d\n', nStrains);
fprintf('Total genes:       %d\n', nGenes);
fprintf('Core genes:        %d (%.1f%%)\n', nCore,      100*nCore/nGenes);
fprintf('Accessory genes:   %d (%.1f%%)\n', nAccessory, 100*nAccessory/nGenes);

% Color map: green (low protection) to blue (high protection)
% Red (low protection) → white (mid) → blue (high protection)
protCmap = [zeros(256,1), linspace(0.8,0,256)', linspace(0,1,256)'];
minP = min(protection);
maxP = max(protection);

%% Figure 1 — Pan-genome pie chart
figure('Position', [100 100 400 400]);
pie([nCore nAccessory]);
colormap([0 0.45 0.74; 0.9 0.5 0.3]);
legend({'Core genome', 'Accessory genome'}, 'Location', 'southoutside', 'FontSize', 14);
title(sprintf('Pan-genome: %d CDS', nGenes), 'FontSize', 16);

%% Figure 2 — PCA of accessory gene matrix colored by protection score
accessoryData          = data(:, accessoryMask);
accessoryData_centered = accessoryData - mean(accessoryData, 1);
[~, score, ~, ~, explained] = pca(accessoryData_centered);

figure('Position', [100 100 600 500]);
hold on;

for i = 1:nStrains
    colorIdx = round((protection(i) - minP) / (maxP - minP) * 255) + 1;
    scatter(score(i,1), score(i,2), 180, protCmap(colorIdx,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

for i = 1:nStrains
    text(score(i,1)+0.3, score(i,2)+0.3, ...
    strrep(strainNames{i}, 'ST1-', 'ST1.'), ...
    'FontSize', 18, 'Interpreter', 'none');
end

colormap(protCmap);
cb = colorbar;
cb.Label.String = 'Protection score';
clim([minP maxP]);
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('PCA: Accessory gene presence/absence');
set(gca, 'FontName', 'Arial', 'FontSize', 20);
box on; hold off;

%% Figure 3 — clcA_1 presence across strains ordered by protection score
clca1_idx = find(strcmp(geneNames, 'clcA_1'));

if isempty(clca1_idx)
    fprintf('WARNING: clcA_1 not found\n');
else
    clca1_sorted = data_sorted(:, clca1_idx);
    clca1_all    = data(:, clca1_idx);
    prot_with    = protection(clca1_all == 1);
    prot_without = protection(clca1_all == 0);
    p_clca1      = ranksum(prot_with, prot_without);

    fprintf('\n=== clcA_1 analysis ===\n');
    fprintf('Strains with clcA_1 (%d): %s\n', sum(clca1_all), ...
        strjoin(strainNames(clca1_all==1), ', '));
    fprintf('Mean protection with:    %.2f\n', mean(prot_with));
    fprintf('Mean protection without: %.2f\n', mean(prot_without));
    fprintf('Wilcoxon p = %.3f\n', p_clca1);

    figure('Position', [100 100 800 280]);
    hold on;

    for i = 1:nStrains
        if clca1_sorted(i) == 1
            scatter(i, 1, 200, [0 0.45 0.74], 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        else
            scatter(i, 1, 200, [1 1 1], 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
    end

    scatter(nan, nan, 100, [0 0.45 0.74], 'filled', 'DisplayName', 'Present');
    scatter(nan, nan, 100, [1 1 1], 'filled', ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Absent');

    set(gca, 'XTick', 1:nStrains);
    set(gca, 'XTickLabel', strrep(strains_sorted, 'ST1-', 'ST1.'));
    set(gca, 'XTickLabelRotation', 45);
    set(gca, 'TickLabelInterpreter', 'none');
    set(gca, 'YTick', []);
    set(gca, 'FontName', 'Arial', 'FontSize', 14);
    xlabel('Strains (ordered by protection score, low to high)');
    title(sprintf('clcA_1 presence across ST1 strains (Wilcoxon p = %.3f)', p_clca1), ...
        'Interpreter', 'none');
    legend('Location', 'northeast');
    xlim([0 nStrains+1]); ylim([0.5 1.5]);
    box on; hold off;
end

%% Compute gene correlations with protection score
fprintf('\nComputing gene correlations — this may take a moment...\n');
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

% Remove invariant genes
validMask         = correlations ~= 0;
validCorr         = correlations(validMask);
validPvals        = pvals_gene(validMask);
validNames        = geneNames(validMask);

% Sort by correlation descending
[validCorr_sorted, wSortIdx] = sort(validCorr, 'descend');
validPvals_sorted = validPvals(wSortIdx);
validNames_sorted = validNames(wSortIdx);
nValid            = length(validCorr_sorted);

% Top 10 and bottom 10
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

clca1_in_plot = find(strcmp(plot_names, 'clcA_1'));

figure('Position', [100 100 950 520]);
hold on;

for i = 1:length(plot_idx)
    bar(i, plot_corr(i), 'FaceColor', barColors(i,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.85);
end

% Highlight clcA_1 with red edge
if ~isempty(clca1_in_plot)
    bar(clca1_in_plot, plot_corr(clca1_in_plot), ...
        'FaceColor', barColors(clca1_in_plot,:), ...
        'EdgeColor', 'r', 'LineWidth', 2.5);
    text(clca1_in_plot, plot_corr(clca1_in_plot) + 0.03, 'clcA\_1', ...
        'FontSize', 11, 'HorizontalAlignment', 'center', ...
        'Color', 'r', 'FontWeight', 'bold', 'Interpreter', 'tex');
end

% Bonferroni threshold lines
yline(r_thresh,  '--k', 'LineWidth', 1.5, ...
    'Label', sprintf('Bonferroni p < %.1e', bonferroni_threshold), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
yline(-r_thresh, '--k', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 0.5);

% Gap between top and bottom 10
xline(10.5, 'k:', 'LineWidth', 1.5);

cleanNames = strrep(plot_names, '_', '\_');
set(gca, 'XTick', 1:length(plot_idx));
set(gca, 'XTickLabel', cleanNames);
set(gca, 'XTickLabelRotation', 45);
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'FontName', 'Arial', 'FontSize', 14);
xlabel('Gene (top 10 positive | bottom 10 negative correlation with protection)');
ylabel('Spearman r with protection score');
title('Gene correlations with protection score');
xlim([0 length(plot_idx)+1]);
box on; grid on; hold off;

fprintf('\n=== Waterfall plot summary ===\n');
fprintf('Bonferroni threshold: p < %.2e (|r| > %.3f)\n', bonferroni_threshold, r_thresh);
fprintf('\nTop 10 positively correlated genes:\n');
for i = 1:10
    sig = '';
    if validPvals_sorted(i) < bonferroni_threshold
        sig = ' *** survives Bonferroni';
    end
    fprintf('  %s: r = %.3f, p = %.4f%s\n', validNames_sorted{i}, ...
        validCorr_sorted(i), validPvals_sorted(i), sig);
end
fprintf('\nBottom 10 negatively correlated genes:\n');
for i = nValid-9:nValid
    sig = '';
    if validPvals_sorted(i) < bonferroni_threshold
        sig = ' *** survives Bonferroni';
    end
    fprintf('  %s: r = %.3f, p = %.4f%s\n', validNames_sorted{i}, ...
        validCorr_sorted(i), validPvals_sorted(i), sig);
end
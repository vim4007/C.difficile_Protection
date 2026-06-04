%% Load data
% tbl = readtable('/Users/vmishra/C.difficile_Protection/data/16s_sequencing/tblAbund.xls');
base_dir = fileparts(fileparts(fileparts(which('microbiome_analysis'))));
tbl = readtable(fullfile(base_dir, 'data', '16s_sequencing', 'tblAbund.xls'));
tbl.Initial_infection   = string(tbl.Initial_infection);
tbl.Initial_antibiotics = string(tbl.Initial_antibiotics);

genusCols = tbl.Properties.VariableNames(6:end);

meanAbund = mean(tbl{:, genusCols}, 1);
[~, sortIdx] = sort(meanAbund, 'descend');
% Exclude pre-existing 'Other' col, take top 9 named taxa, combine rest + Other into one
namedIdx = sortIdx(~strcmp(genusCols(sortIdx), 'Other'));
top9     = genusCols(namedIdx(1:9));
otherIdx = namedIdx(10:end);
tbl.Other_combined = tbl.Other + sum(tbl{:, genusCols(otherIdx)}, 2);
plotCols   = [top9, {'Other_combined'}];
plotLabels = [top9, {'Other'}];
nPlot      = length(plotCols);

groupKeys   = {'uninfected', 'ST1_75',    'ST1_12'};
groupLabels = {'Uninfected', 'ST1.75',    'ST1.12'};
groupColors = [0.5 0.5 0.5; 0 0.45 0.74; 0.9 0.5 0.3];
analysisDays = [0 1 2 7];

hexColors = {
    '#5875DE', ... %  1 Lactobacillus              
    '#1CF8EC', ... %  2 Bacteroides                
    '#2A9D8F', ... %  3 unclassified_Muribaculaceae_ASV_4 
    '#0D7E2B', ... %  4 Enterococcus               
    '#CA0BE8', ... %  5 Akkermansia                
    '#FBA22E', ... %  6 Turicibacter               
    '#BEA89A', ... %  7 LachnospiraceaeNK4A136Group 
    '#7D6E65', ... %  8 Clostridioides             
    '#1A7A73', ... %  9 unclassified_Muribaculaceae_ASV_6 
    '#D4B896', ... % 10 Other_combined             
};
genusPalette = reshape(cell2mat(cellfun(@(h) sscanf(h(2:end),'%2x%2x%2x')'/255, hexColors, 'UniformOutput', false)), 3, [])';

shannon = @(p) -sum(p(p>0) .* log(p(p>0)));
tbl.Shannon = zeros(height(tbl), 1);
for i = 1:height(tbl)
    p = tbl{i, genusCols} / 100;
    tbl.Shannon(i) = shannon(p);
end

%% Figure 1 — Stacked bar Day 0: neg vs mnvc averaged across all mice
day0      = tbl(tbl.Day == 0, :);
neg_sub   = day0(day0.Initial_antibiotics == 'neg',  :);
mnvc_sub  = day0(day0.Initial_antibiotics == 'mnvc', :);
neg_mean  = mean(neg_sub{:,  plotCols}, 1) / 100;
mnvc_mean = mean(mnvc_sub{:, plotCols}, 1) / 100;

figure('Position', [100 100 420 500]);
b = bar([neg_mean; mnvc_mean], 'stacked', 'EdgeColor', 'none');
for i = 1:nPlot; b(i).FaceColor = genusPalette(i,:); end
set(gca, 'XTick', 1:2, ...
    'XTickLabel', {'No antibiotic', 'Antibiotic (mnvc)'}, ...
    'XTickLabelRotation', 60, 'FontSize', 40);
ylabel('Relative abundance');
title('Day 0 — Antibiotic effect');
ylim([0 1]);
lgd = legend(plotLabels, 'Location', 'eastoutside');
set(lgd, 'Interpreter', 'none', 'FontSize', 30);
box on;

%% Figures 2-4 — Stacked bars Days 1, 2, 7 — mnvc only averaged per group
postDays  = [1 2 7];
dayTitles = {'Day 1', 'Day 2', 'Day 7'};

for d = 1:length(postDays)
    sub = tbl(tbl.Day == postDays(d) & tbl.Initial_antibiotics == 'mnvc', :);
    barData = zeros(3, nPlot);
    for g = 1:3
        grp = sub(sub.Initial_infection == groupKeys{g}, :);
        if ~isempty(grp)
            barData(g,:) = mean(grp{:, plotCols}, 1) / 100;
        end
    end
    figure('Position', [100 100 450 500]);
    b = bar(barData, 'stacked', 'EdgeColor', 'none');
    for i = 1:nPlot; b(i).FaceColor = genusPalette(i,:); end
    set(gca, 'XTick', 1:3, 'XTickLabel', groupLabels, ...
        'XTickLabelRotation', 30, 'FontSize', 40);
    ylabel('Relative abundance');
    title(dayTitles{d});
    ylim([0 1]);
    lgd = legend(plotLabels, 'Location', 'eastoutside');
    set(lgd, 'Interpreter', 'none', 'FontSize', 30);
    box on;
end

%% Figure 5 — Alpha diversity (Shannon) mnvc only ST1.75 vs ST1.12
figure('Position', [100 100 800 500]);
hold on;

barW = 0.15;
xPos = 1:length(analysisDays);

for d = 1:length(analysisDays)
    day = analysisDays(d);
    sub = tbl(tbl.Day == day & tbl.Initial_antibiotics == 'mnvc', :);

    d75 = sub(sub.Initial_infection == 'ST1_75', :);
    d12 = sub(sub.Initial_infection == 'ST1_12', :);

    for g = 1:2
        if g == 1; grp = d75; else; grp = d12; end
        if isempty(grp), continue; end
        x = xPos(d) + (g-1.5)*barW*2;
        scatter(repmat(x, height(grp), 1), grp.Shannon, 200, ...
            groupColors(g+1,:), 'filled', 'MarkerFaceAlpha', 0.8, ...
            'HandleVisibility', 'off');
        errorbar(x, mean(grp.Shannon), std(grp.Shannon), 'k.', ...
            'LineWidth', 2.5, 'CapSize', 8, 'HandleVisibility', 'off');
    end

    if height(d75) >= 2 && height(d12) >= 2
        try
            p = ranksum(d75.Shannon, d12.Shannon);
            text(xPos(d), max([d75.Shannon; d12.Shannon]) + 0.05, ...
                sprintf('p=%.2f', p), 'HorizontalAlignment', 'center', 'FontSize', 20);
        catch
        end
    end
end

scatter(nan, nan, 200, groupColors(2,:), 'filled', 'DisplayName', 'ST1.75');
scatter(nan, nan, 200, groupColors(3,:), 'filled', 'DisplayName', 'ST1.12');

set(gca, 'XTick', xPos, ...
    'XTickLabel', {'Day 0', 'Day 1', 'Day 2', 'Day 7'}, 'FontSize', 30);
xlabel('Day');
ylabel('Shannon diversity index');
%title('Alpha diversity');
legend('Location', 'northeast');
box on; grid on; hold off;

%% Figures 6-9 — PCoA one per day
rng(42);

for d = 1:length(analysisDays)
    day = analysisDays(d);

    if day == 0
     
        sub = tbl(tbl.Day == day, :);
        abd = sub{:, genusCols} / 100;
        nS  = height(sub);

        BC = zeros(nS);
        for i = 1:nS
            for j = i+1:nS
                xi = abd(i,:); xj = abd(j,:);
                v  = sum(abs(xi-xj)) / (sum(xi)+sum(xj));
                BC(i,j) = v; BC(j,i) = v;
            end
        end

        [coords, eigvals] = cmdscale(BC);
        varExp = eigvals / sum(eigvals(eigvals>0)) * 100;

        figure('Position', [100 100 500 450]);
        hold on;

        neg_idx  = sub.Initial_antibiotics == 'neg';
        mnvc_idx = sub.Initial_antibiotics == 'mnvc';

        scatter(coords(neg_idx,1),  coords(neg_idx,2),  180, ...
            [0.85 0.33 0.10], 'filled', 'MarkerFaceAlpha', 0.8, ...
            'DisplayName', 'No antibiotic');
        scatter(coords(mnvc_idx,1), coords(mnvc_idx,2), 180, ...
            [0 0.45 0.74], 'filled', 'MarkerFaceAlpha', 0.8, ...
            'DisplayName', 'Antibiotic (mnvc)');

        xlabel(sprintf('PCo1 (%.1f%%)', varExp(1)));
        ylabel(sprintf('PCo2 (%.1f%%)', varExp(2)));
        title('PCoA Bray-Curtis — Day 0');
        legend('Location', 'bestoutside');
        set(gca, 'FontSize', 20); box on; hold off;

    else
        sub = tbl(tbl.Day == day & tbl.Initial_antibiotics == 'mnvc' & ...
            (tbl.Initial_infection == 'ST1_75' | ...
             tbl.Initial_infection == 'ST1_12'), :);
        abd = sub{:, genusCols} / 100;
        nS  = height(sub);

        BC = zeros(nS);
        for i = 1:nS
            for j = i+1:nS
                xi = abd(i,:); xj = abd(j,:);
                v  = sum(abs(xi-xj)) / (sum(xi)+sum(xj));
                BC(i,j) = v; BC(j,i) = v;
            end
        end

        [coords, eigvals] = cmdscale(BC);
        varExp = eigvals / sum(eigvals(eigvals>0)) * 100;

        figure('Position', [100 100 500 450]);
        hold on;

        idx75 = sub.Initial_infection == 'ST1_75';
        idx12 = sub.Initial_infection == 'ST1_12';

        scatter(coords(idx75,1), coords(idx75,2), 180, ...
            groupColors(2,:), 'filled', 'MarkerFaceAlpha', 0.8, ...
            'DisplayName', 'ST1.75');
        scatter(coords(idx12,1), coords(idx12,2), 180, ...
            groupColors(3,:), 'filled', 'MarkerFaceAlpha', 0.8, ...
            'DisplayName', 'ST1.12');

        xlabel(sprintf('PCo1 (%.1f%%)', varExp(1)));
        ylabel(sprintf('PCo2 (%.1f%%)', varExp(2)));
        title(sprintf('PCoA Bray-Curtis — Day %d', day));
        legend('Location', 'bestoutside');
        set(gca, 'FontSize', 20); box on; hold off;
    end
end

%% Figure 10 — Clostridioides abundance mnvc only ST1.75 vs ST1.12
figure('Position', [100 100 700 500]);
hold on;
barW = 0.15;

for d = 1:length(analysisDays)
    day = analysisDays(d);
    sub = tbl(tbl.Day == day & tbl.Initial_antibiotics == 'mnvc', :);

    d75 = sub(sub.Initial_infection == 'ST1_75', :);
    d12 = sub(sub.Initial_infection == 'ST1_12', :);

    for g = 1:2
        if g == 1; grp = d75; else; grp = d12; end
        if isempty(grp), continue; end
        x = d + (g-1.5)*barW*2;
        scatter(repmat(x, height(grp), 1), grp.Clostridioides, 200, ...
            groupColors(g+1,:), 'filled', 'MarkerFaceAlpha', 0.7, ...
            'HandleVisibility', 'off');
        errorbar(x, mean(grp.Clostridioides,'omitnan'), ...
            std(grp.Clostridioides,'omitnan'), 'k.', ...
            'LineWidth', 1.5, 'CapSize', 5, 'HandleVisibility', 'off');
    end

    if height(d75) >= 2 && height(d12) >= 2
        try
            p = ranksum(d75.Clostridioides, d12.Clostridioides);
            fprintf('  Day %d: ST1.75=%.2f%%, ST1.12=%.2f%%, p=%.3f\n', ...
                day, mean(d75.Clostridioides,'omitnan'), ...
                mean(d12.Clostridioides,'omitnan'), p);
        catch
        end
    end
end

% Legend
scatter(nan, nan, 200, groupColors(2,:), 'filled', 'DisplayName', 'ST1.75');
scatter(nan, nan, 200, groupColors(3,:), 'filled', 'DisplayName', 'ST1.12');

set(gca, 'XTick', 1:length(analysisDays), ...
    'XTickLabel', {'Day 0', 'Day 1', 'Day 2', 'Day 7'}, 'FontSize', 25);
xlabel('Day');
ylabel('Clostridioides relative abundance (%)');
title('Clostridioides abundance (mnvc only)');
legend('Location', 'northeast');
box on; grid on; hold off;
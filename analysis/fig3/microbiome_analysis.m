
%% Load data
tbl = readtable('/Users/vmishra/C.difficile_Protection/data/16s_sequencing/tblAbund.xls');
tbl.Initial_infection   = string(tbl.Initial_infection);
tbl.Initial_antibiotics = string(tbl.Initial_antibiotics);

genusCols = tbl.Properties.VariableNames(6:end);

meanAbund = mean(tbl{:, genusCols}, 1);
[~, sortIdx] = sort(meanAbund, 'descend');
top12    = genusCols(sortIdx(1:12));
otherIdx = sortIdx(13:end);
tbl.Other_combined = sum(tbl{:, genusCols(otherIdx)}, 2);
plotCols   = [top12, {'Other_combined'}];
plotLabels = [top12, {'Other'}];
nPlot      = length(plotCols);

groupKeys   = {'uninfected', 'ST1_75',    'ST1_12'};
groupLabels = {'Uninfected', 'ST1.75',    'ST1.12'};
groupColors = [0.5 0.5 0.5; 0 0.45 0.74; 0.9 0.5 0.3];
analysisDays = [0 1 2 7];
genusPalette = turbo(nPlot);

shannon = @(p) -sum(p(p>0) .* log(p(p>0)));
tbl.Shannon = zeros(height(tbl), 1);
for i = 1:height(tbl)
    p = tbl{i, genusCols} / 100;
    tbl.Shannon(i) = shannon(p);
end

%% 
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
    'XTickLabelRotation', 60, 'FontSize', 30);
ylabel('Relative abundance');
title('Day 0 — Antibiotic effect');
ylim([0 1]);
lgd = legend(plotLabels, 'Location', 'eastoutside');
set(lgd, 'Interpreter', 'none', 'FontSize', 20);
box on;

%% 
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
        'XTickLabelRotation', 30, 'FontSize', 30);
    ylabel('Relative abundance');
    title(dayTitles{d});
    ylim([0 1]);
    lgd = legend(plotLabels, 'Location', 'eastoutside');
    set(lgd, 'Interpreter', 'none', 'FontSize', 20);
    box on;
end

%% 
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
        scatter(repmat(x, height(grp), 1), grp.Shannon, 150, ...
            groupColors(g+1,:), 'filled', 'MarkerFaceAlpha', 0.7, ...
            'HandleVisibility', 'off');
        errorbar(x, mean(grp.Shannon), std(grp.Shannon), 'k.', ...
            'LineWidth', 1.5, 'CapSize', 5, 'HandleVisibility', 'off');
    end

    if ~isempty(d75) && ~isempty(d12)
        p = ranksum(d75.Shannon, d12.Shannon);
       
    end
end

scatter(nan, nan, 150, groupColors(2,:), 'filled', 'DisplayName', 'ST1.75');
scatter(nan, nan, 150, groupColors(3,:), 'filled', 'DisplayName', 'ST1.12');

set(gca, 'XTick', xPos, ...
    'XTickLabel', {'Day 0', 'Day 1', 'Day 2', 'Day 7'}, 'FontSize', 30);
xlabel('Day');
ylabel('Shannon diversity index');
title('Alpha diversity (mnvc only)');
legend('Location', 'northeast');
box on; grid on; hold off;

%% 
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

%% 
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
        scatter(repmat(x, height(grp), 1), grp.Clostridioides, 40, ...
            groupColors(g+1,:), 'filled', 'MarkerFaceAlpha', 0.7, ...
            'HandleVisibility', 'off');
        errorbar(x, mean(grp.Clostridioides,'omitnan'), ...
            std(grp.Clostridioides,'omitnan'), 'k.', ...
            'LineWidth', 1.5, 'CapSize', 5, 'HandleVisibility', 'off');
    end

    if ~isempty(d75) && ~isempty(d12)
        p = ranksum(d75.Clostridioides, d12.Clostridioides);
        fprintf('  Day %d: ST1.75=%.2f%%, ST1.12=%.2f%%, p=%.3f\n', ...
            day, mean(d75.Clostridioides,'omitnan'), ...
            mean(d12.Clostridioides,'omitnan'), p);
    end
end

% Legend
scatter(nan, nan, 40, groupColors(2,:), 'filled', 'DisplayName', 'ST1.75');
scatter(nan, nan, 40, groupColors(3,:), 'filled', 'DisplayName', 'ST1.12');

set(gca, 'XTick', 1:length(analysisDays), ...
    'XTickLabel', {'Day 0', 'Day 1', 'Day 2', 'Day 7'}, 'FontSize', 20);
xlabel('Day');
ylabel('Clostridioides relative abundance (%)');
title('Clostridioides abundance (mnvc only)');
legend('Location', 'northeast');
box on; grid on; hold off;

%% --- Load data ---
% B   = readtable('/Users/vmishra/C.difficile_Protection/data/Biolog/Biolog_growth_matrix.xlsx','VariableNamingRule','preserve');
% G   = readtable('/Users/vmishra/C.difficile_Protection/data/Biolog/strain_groups.xlsx',           'VariableNamingRule','preserve');
% MoA = readtable('/Users/vmishra/C.difficile_Protection/data/Biolog/moas.xlsx',                    'VariableNamingRule','preserve');

base_dir = fileparts(fileparts(mfilename('fullpath')));
B   = readtable(fullfile(base_dir, 'data', 'Biolog', 'Biolog_growth_matrix.xlsx'), 'VariableNamingRule', 'preserve');
G   = readtable(fullfile(base_dir, 'data', 'Biolog', 'strain_groups.xlsx'),        'VariableNamingRule', 'preserve');
MoA = readtable(fullfile(base_dir, 'data', 'Biolog', 'moas.xlsx'),                 'VariableNamingRule', 'preserve');

B = B(~strcmp(B.Metabolites,'Negative Control'),:);
G.Strains_norm = strrep(G.Strains,'-','_');
G_ST1 = G(startsWith(G.Strains_norm,'ST1'),:);
st1   = G_ST1.Strains_norm;
M     = table2array(B(:,st1));         
prot  = G_ST1.Protection_Estimate;
n_met = height(B);  n_st1 = numel(st1);

MoA.category = regexprep(MoA.MoA,'^C-Source, ','');
[~, loc] = ismember(B.Metabolites, MoA.Chemical);
cats_per_row = MoA.category(loc);

cmap_gb = [linspace(0.20,0.10,256)', linspace(0.70,0.30,256)', linspace(0.30,0.75,256)'];

%% 
p_w   = nan(n_met,1);
rho_m = nan(n_met,1);
n_g   = sum(M,2);

for i = 1:n_met
    g = M(i,:)';
    if sum(g)<2 || sum(~g)<2, continue; end
    p_w(i)   = ranksum(prot(g==1), prot(g==0));
    rho_m(i) = corr(double(g), prot, 'Type','Spearman');
end

m_tests  = sum(~isnan(p_w));
bonf_thr = 0.05 / m_tests;
sig_bonf = p_w < bonf_thr;

W = table(B.Metabolites, n_g, rho_m, p_w, sig_bonf, cats_per_row, ...
    'VariableNames',{'Metabolite','n_growers','rho','p_wilcox','sig_bonf','category'});
W = sortrows(W,'p_wilcox');

%% 
cats = unique(cats_per_row(~cellfun(@isempty,cats_per_row)));
rho = nan(numel(cats),1); p_c = nan(numel(cats),1); tot = nan(numel(cats),1);
for j = 1:numel(cats)
    r = strcmp(cats_per_row, cats{j});
    tot(j) = sum(r);
    counts = sum(M(r,:),1)';
    if std(counts)>0
        [rho(j), p_c(j)] = corr(counts, prot, 'Type','Spearman');
    end
end

m_cat        = sum(~isnan(p_c));
bonf_thr_cat = 0.05 / m_cat;
sig_bonf_c   = p_c < bonf_thr_cat;

MR = table(cats, tot, rho, p_c, sig_bonf_c, ...
    'VariableNames',{'Category','TotalMets','rho','p','sig_bonf'});
MR = sortrows(MR,'rho','descend');

%% 
cat_pal = lines(numel(cats));
cat_map = containers.Map(cats, num2cell(cat_pal,2));

MR2 = sortrows(MR,'rho','descend');   % highest rho on the left
figure('Position',[100 100 750 500]);

for i = 1:height(MR2)
    c = cat_map(MR2.Category{i});
    bar(i, MR2.rho(i),'FaceColor',c,'EdgeColor','k'); hold on;

    if MR2.p(i) < bonf_thr_cat
        mk = '**';
    elseif MR2.p(i) < 0.05
        mk = '';
    else
        mk = '';
    end
    if ~isempty(mk)
        y_off = 0.02 * sign(MR2.rho(i));
        va = 'bottom'; if MR2.rho(i)<0, va='top'; end
        text(i, MR2.rho(i)+y_off, mk, 'FontSize',14,'FontWeight','bold', ...
             'HorizontalAlignment','center','VerticalAlignment',va);
    end
end

set(gca,'XTick',1:height(MR2), ...
    'XTickLabel',arrayfun(@(i) sprintf('%s (%d)', MR2.Category{i}, MR2.TotalMets(i)), ...
                          1:height(MR2),'UniformOutput',false), ...
    'XTickLabelRotation',45,'FontSize',20);
yline(0,'k');
ylabel('Spearman \rho vs protection');
ylim([0 1]);
title('Chemical class associations with protection');
box off;

h = arrayfun(@(j) patch(NaN,NaN,cat_map(MR2.Category{j})), 1:height(MR2));
legend(h, MR2.Category,'Location','northeastoutside','Box','off');

%% 
[~, ord] = sort(prot,'descend');
prot_ord = prot(ord);
pn = (prot_ord - min(prot)) / (max(prot)-min(prot));
strain_colors = cmap_gb(max(1,round(pn*255)+1),:);

cat_counts = zeros(numel(cats),1);
for j = 1:numel(cats)
    cat_counts(j) = sum(strcmp(cats_per_row, cats{j}));
end

figure('Position',[100 100 900 550]);
for j = 1:numel(cats)
    r = strcmp(cats_per_row, cats{j});
    col_counts = sum(M(r,:),1);
    col_counts = col_counts(ord);
    total = sum(col_counts);
    if total == 0, continue; end
    props = col_counts / total;
    y0 = 0;
    for k = 1:n_st1
        rectangle('Position',[j-0.4, y0, 0.8, props(k)], ...
            'FaceColor',strain_colors(k,:),'EdgeColor','w','LineWidth',0.3); hold on;
        y0 = y0 + props(k);
    end
    text(j, 1.03, sprintf('%d', cat_counts(j)), ...
        'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
end
xtick_labels = arrayfun(@(j) sprintf('%s (%d)', cats{j}, cat_counts(j)), ...
                        1:numel(cats),'UniformOutput',false);

set(gca,'XTick',1:numel(cats),'XTickLabel',xtick_labels, ...
    'XTickLabelRotation',45,'FontSize',14);
ylabel('Proportion'); ylim([0 1.10]); xlim([0.5 numel(cats)+0.5]);
title('Stacked strain contributions per chemical class');
colormap(cmap_gb); cb = colorbar; cb.Label.String='Protection estimate';
caxis([min(prot) max(prot)]); box off;

%% 
X = double(M');
[coeff, scores, ~, ~, ev] = pca(X);
fprintf('\nPCA variance: PC1 %.1f%%, PC2 %.1f%%\n', ev(1), ev(2));

figure('Position',[100 100 650 520]);
colormap(cmap_gb);
scatter(scores(:,1),scores(:,2),250,prot,'filled','MarkerEdgeColor','k'); hold on;
for i = 1:n_st1
    text(scores(i,1)+0.1, scores(i,2)+0.08, strrep(st1{i},'ST1_','ST1-'),'FontSize',20);
end
if any(strcmp(B.Properties.VariableNames,'VPI'))
    vpi = (double(table2array(B(:,'VPI'))') - mean(X,1)) * coeff;
    scatter(vpi(1),vpi(2),250,'r^','filled','MarkerEdgeColor','k');
    text(vpi(1)+0.15, vpi(2), 'VPI','Color','r','FontWeight','bold','FontSize',20);
end
xline(0,'--','Color',[.6 .6 .6]); yline(0,'--','Color',[.6 .6 .6]);
xlabel(sprintf('PC1 (%.1f%%)',ev(1))); ylabel(sprintf('PC2 (%.1f%%)',ev(2)));
%title('PCA of BIOLOG growth profiles');
cb = colorbar; cb.Label.String='Protection estimate';
clim([min(prot) max(prot)]); box off;
set(gca,'FontSize',30);

%%

moa_counts = zeros(n_st1, numel(cats));
for j = 1:numel(cats)
    r = strcmp(cats_per_row, cats{j});
    moa_counts(:,j) = sum(M(r,:),1)';
end

[prot_sorted, ord] = sort(prot,'descend');
moa_sorted = moa_counts(ord,:);
st1_sorted = st1(ord);

figure('Position',[100 100 1100 550]);

yyaxis left
b = bar(1:n_st1, moa_sorted, 'stacked', 'EdgeColor','w','LineWidth',0.5);
for j = 1:numel(cats)
    b(j).FaceColor = cat_map(cats{j});
end
ylabel('Number of carbon sources supporting growth');
set(gca,'YColor','k');
yyaxis right
x_vals = (1:n_st1)';

plot(x_vals, prot_sorted, 'k-o', 'LineWidth',1.5, ...
    'MarkerSize',6,'MarkerFaceColor','k');
ylabel('Protection estimate');
set(gca,'YColor','k');

hold on;
p_fit = polyfit(x_vals, prot_sorted, 1);
y_fit = polyval(p_fit, x_vals);
plot(x_vals, y_fit, 'r--', 'LineWidth',1.8);
niche_breadth = sum(moa_sorted, 2);
[rho, pval] = corr(niche_breadth, prot_sorted, 'Type','Spearman');

if pval < 0.001
    p_str = 'p < 0.001';
else
    p_str = sprintf('p = %.3f', pval);
end

ax = gca;
y_range = ax.YLim;
annotation_str = sprintf('\\rho = %.2f, %s', rho, p_str);
text(n_st1 * 0.97, y_range(1) + 0.90*(y_range(2)-y_range(1)), ...
    annotation_str, 'FontSize',20,'HorizontalAlignment','right', ...
    'Color','r','FontWeight','bold');

set(gca,'XTick',1:n_st1, ...
    'XTickLabel',strrep(st1_sorted,'ST1_','ST1-'), ...
    'XTickLabelRotation',45,'FontSize',20);
xlabel('Strains sorted by protection (high \rightarrow low)');
xlim([0.5 n_st1+0.5]);
cat_counts = arrayfun(@(j) sum(strcmp(cats_per_row, cats{j})), 1:numel(cats));
leg_labels = arrayfun(@(j) sprintf('%s (%d)', cats{j}, cat_counts(j)), ...
                      1:numel(cats),'UniformOutput',false);
legend(b, leg_labels,'Location','northeastoutside','Box','off','FontSize',20);

box off;
%%

grp_labels = G_ST1.Protection_Group; 
grp_num = zeros(n_st1,1);
grp_num(strcmp(grp_labels,'High'))   = 3;
grp_num(strcmp(grp_labels,'Medium')) = 2;
grp_num(strcmp(grp_labels,'Low'))    = 1;

X = M';  

n_correct = 0;
rf_importance = zeros(1, n_met);
for i = 1:n_st1
    train_idx = setdiff(1:n_st1, i);
    mdl = TreeBagger(200, X(train_idx,:), grp_num(train_idx), ...
        'Method','classification','OOBPredictorImportance','off');
    pred = str2double(predict(mdl, X(i,:)));
    if pred == grp_num(i); n_correct = n_correct + 1; end

    full_mdl = TreeBagger(200, X, grp_num, 'Method','classification', ...
        'OOBPredictorImportance','on');
    rf_importance = rf_importance + full_mdl.OOBPermutedPredictorDeltaError;
end
rf_importance = rf_importance / n_st1;
accuracy = n_correct / n_st1 * 100;
fprintf('RF accuracy: %.0f%%\n', accuracy);

[~, rf_ord] = sort(rf_importance,'descend');
fprintf('Top 5 RF features:\n');
for i = 1:5
    fprintf('  %s\n', B.Metabolites{rf_ord(i)});
end

% 
[B_lasso, fitinfo] = lasso(X, prot, 'Alpha',1,'NumLambda',10,'CV',5);
idx_min = fitinfo.IndexMinMSE;
coef = B_lasso(:, idx_min);
selected = find(coef ~= 0);
fprintf('\nLasso selected metabolites:\n');
for i = 1:numel(selected)
    fprintf('  %s  (coef = %.3f)\n', B.Metabolites{selected(i)}, coef(selected(i)));
end

%%
gly_rank_rf = find(rf_ord == find(strcmp(B.Metabolites,'Glyoxylic acid')));
fprintf('\nGlyoxylic acid RF rank: %d of %d\n', gly_rank_rf, n_met);
fprintf('Glyoxylic acid Lasso coef: %.3f\n', coef(strcmp(B.Metabolites,'Glyoxylic acid')));

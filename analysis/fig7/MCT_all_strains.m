
% tbl = readtable('Biolog_growth_matrix.xlsx', 'VariableNamingRule', 'preserve');

script_dir = fileparts(which('MCT_all_strains'));
base_dir = fileparts(fileparts(script_dir));
tbl = readtable(fullfile(base_dir, 'data', 'Biolog', 'Biolog_growth_matrix.xlsx'), 'VariableNamingRule', 'preserve');

metabolites  = tbl{:, 1};
strain_names = tbl.Properties.VariableNames(2:end);
data         = tbl{:, 2:end};

neg_idx = strcmp(metabolites, 'Negative Control');
data(neg_idx, :)     = [];
metabolites(neg_idx) = [];

vpi_col = find(strcmp(strain_names, 'VPI'));
if isempty(vpi_col)
    error('VPI column not found in Biolog file.');
end
vpi_vec   = data(:, vpi_col);
st1_cols  = setdiff(1:numel(strain_names), vpi_col);
st1_names = strain_names(st1_cols);   % e.g. 'ST1_75'
n_strains = numel(st1_cols);

%%
prot_tbl = readtable('strain_groups.xlsx', 'VariableNamingRule', 'preserve');

prot_names  = strrep(prot_tbl.Strains, '-', '_');
prot_scores = prot_tbl.Protection_Estimate;
prot_groups = prot_tbl.Protection_Group;  

%% 
n1 = zeros(n_strains, 1);
n2 = zeros(n_strains, 1);
ns = zeros(n_strains, 1);

for i = 1:n_strains
    st1_vec = data(:, st1_cols(i));
    n1(i) = sum(st1_vec == 1 & vpi_vec == 0);
    n2(i) = sum(st1_vec == 0 & vpi_vec == 1);
    ns(i) = sum(st1_vec == 1 & vpi_vec == 1);
end

%% 

p.U      = 0.0;

kappa = (n1 - n2) + p.U .* ns;
S     = 2 .* n1 .* n2 ./ (n1 + n2);

outcome = strings(n_strains, 1);
for i = 1:n_strains
    if     kappa(i) >  S(i);  outcome(i) = "ST1 excludes VPI";
    elseif kappa(i) < -S(i);  outcome(i) = "VPI excludes ST1";
    else;                      outcome(i) = "Coexistence";
    end
end



%% 
col_high   = [0.15 0.35 0.75];  
col_medium = [0.60 0.60 0.60];   
col_low    = [0.20 0.65 0.25];   
col_none   = [0.85 0.85 0.85];   

strain_colors = repmat(col_none, n_strains, 1);
strain_group  = repmat("N/A", n_strains, 1);

for i = 1:n_strains
    idx = find(strcmp(prot_names, st1_names{i}));
    if ~isempty(idx)
        strain_group(i) = string(prot_groups{idx});
        switch prot_groups{idx}
            case 'High';   strain_colors(i,:) = col_high;    
            case 'Medium'; strain_colors(i,:) = col_medium;
            case 'Low';    strain_colors(i,:) = col_low;     
        end
    end
end

%%
fig = figure('Units','inches','Position',[1 1 9.5 7],'Color','w');
hold on;

S_max   = max(S) * 1.3;
S_range = linspace(0, S_max, 300);


fill([0 S_range S_max 0],           [0  S_range  0 0],              [0.88 0.93 0.88],'EdgeColor','none','FaceAlpha',0.50);
fill([0 S_range S_max 0],           [0 -S_range  0 0],              [0.93 0.88 0.88],'EdgeColor','none','FaceAlpha',0.50);
fill([0 S_range fliplr(S_range) 0], [0  S_range -fliplr(S_range) 0],[0.87 0.91 0.97],'EdgeColor','none','FaceAlpha',0.60);


plot(S_range,  S_range, 'k--', 'LineWidth', 2.5);
plot(S_range, -S_range, 'k--', 'LineWidth', 2.5);
yline(0,'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',1.0);

for i = 1:n_strains
    scatter(S(i), kappa(i), 180, strain_colors(i,:), 'o', 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth', 0.8);
    lbl = strrep(st1_names{i}, '_', '.');   
    text(S(i) + 0.18, kappa(i), lbl, ...
        'FontSize', 25, 'Color', strain_colors(i,:) * 0.70,'FontWeight','bold');
end


x_lbl = S_max * 0.76;
text(x_lbl,  S_max*0.62, 'ST1 excludes VPI','FontSize',30, ...
    'Color',[0.15 0.42 0.15],'FontWeight','bold','HorizontalAlignment','center');
text(x_lbl, -S_max*0.62, 'VPI excludes ST1','FontSize',30, ...
    'Color',[0.55 0.10 0.10],'FontWeight','bold','HorizontalAlignment','center');
text(S_max*0.20, 0.5, 'Coexistence','FontSize',30, ...
    'Color',[0.15 0.25 0.55],'FontWeight','bold','HorizontalAlignment','center');

h_high = scatter(nan, nan, 80, col_high,   'o','filled','MarkerEdgeColor','k','LineWidth',0.8);
h_med  = scatter(nan, nan, 80, col_medium, 'o','filled','MarkerEdgeColor','k','LineWidth',0.8);
h_low  = scatter(nan, nan, 80, col_low,    'o','filled','MarkerEdgeColor','k','LineWidth',0.8);
legend([h_high h_med h_low], {'High protection','Medium protection','Low protection'}, ...
    'Location','northwest','Box','off','FontSize',20);


xlabel('Stabilizing difference ', 'FontSize',13);
ylabel('Fitness difference ',  'FontSize',13);
%title('MCT Phase Plane: ST1 strains vs VPI10463','FontSize',14,'FontWeight','bold');
xlim([0 S_max]);
ylim([-S_max*1.1  S_max*1.1*2]);
set(gca,'FontSize',35,'Box','off');
grid on;

cmap_discrete = [col_high; col_medium; col_low];
colormap(fig, cmap_discrete);
clim([0 3]);
cb = colorbar('Location','eastoutside');
cb.Ticks      = [0.5 1.5 2.5];
cb.TickLabels = {'High','Medium','Low'};
cb.Label.String   = 'Protection Group';
cb.Label.FontSize = 20;
cb.FontSize       = 20;
cb.TickLength     = 0;


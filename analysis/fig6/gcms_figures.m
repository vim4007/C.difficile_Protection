ext_file = '/Users/vmishra/C.difficile_Protection/data/GCMS/secretome.xlsx';
int_file = '/Users/vmishra/C.difficile_Protection/data/GCMS/intra.xlsx';

ext_table = readtable(ext_file, 'Sheet', 1, 'VariableNamingRule', 'preserve');
int_table = readtable(int_file, 'Sheet', 1, 'VariableNamingRule', 'preserve');

raw_colnames = ext_table.Properties.VariableNames(2:end);
metabolites  = cellfun(@(x) regexprep(x, '_\d+Results$', ''), raw_colnames, 'UniformOutput', false);

sample_ids_ext = ext_table{:,1};
data_ext       = ext_table{:,2:end};

data_ext = data_ext ./ sum(data_ext, 2);

for i = 1:size(data_ext,2)
    col = data_ext(:,i);
    min_nonzero = min(col(col > 0));
    if ~isempty(min_nonzero)
        col(col == 0) = min_nonzero / 2;
    end
    data_ext(:,i) = col;
end

data_ext = log2(data_ext);

valid_ext = std(data_ext, 0, 1) > 0;
data_ext  = data_ext(:, valid_ext);

media_ext      = data_ext(contains(sample_ids_ext,'Media'),  :);
st75_ext       = data_ext(contains(sample_ids_ext,'ST1-75'), :);
st68_ext       = data_ext(contains(sample_ids_ext,'ST1-68'), :);
vpi_ext        = data_ext(contains(sample_ids_ext,'VPI'),    :);

metabolites_ext = metabolites(valid_ext);

sample_ids_int = int_table{:,1};
data_int       = int_table{:,2:end};

data_int = data_int ./ sum(data_int, 2);
for i = 1:size(data_int,2)
    col = data_int(:,i);
    min_nonzero = min(col(col > 0));
    if ~isempty(min_nonzero)
        col(col == 0) = min_nonzero / 2;
    end
    data_int(:,i) = col;
end

data_int = log2(data_int);

valid_int = std(data_int, 0, 1) > 0;
data_int  = data_int(:, valid_int);

media_int      = data_int(contains(sample_ids_int,'Media'),  :);
st75_int       = data_int(contains(sample_ids_int,'ST1-75'), :);
st68_int       = data_int(contains(sample_ids_int,'ST1-68'), :);
vpi_int        = data_int(contains(sample_ids_int,'VPI'),    :);

metabolites_int = metabolites(valid_int);

eps_val = 1e-10; 
alpha   = 0.05;

col75  = [0.20 0.53 0.74];   % blue
col68  = [0.20 0.63 0.29];   % green
colVPI = [0.85 0.15 0.15];   % red

%% 

media_ext_mean = mean(media_ext, 1);
media_int_mean = mean(media_int, 1);

[p_ext_75v68_all, p_ext_75vvpi_all, p_ext_68vvpi_all] = run_ttests(st75_ext, st68_ext, vpi_ext);
[p_int_75v68_all, p_int_75vvpi_all, p_int_68vvpi_all] = run_ttests(st75_int, st68_int, vpi_int);

sig_ext = (p_ext_75v68_all < alpha) | (p_ext_75vvpi_all < alpha) | (p_ext_68vvpi_all < alpha);
sig_int = (p_int_75v68_all < alpha) | (p_int_75vvpi_all < alpha) | (p_int_68vvpi_all < alpha);

fc_75_ext = mean(st75_ext,1) - media_ext_mean;
fc_75_int = mean(st75_int,1) - media_int_mean;

fc_75_ext_sig = fc_75_ext(sig_ext);
fc_75_int_sig = fc_75_int(sig_int);

[~, ext_sort] = sort(fc_75_ext_sig, 'descend');
[~, int_sort] = sort(fc_75_int_sig, 'descend');

ext_idx_sig = find(sig_ext);
int_idx_sig = find(sig_int);

ext_idx = ext_idx_sig(ext_sort);
int_idx = int_idx_sig(int_sort);

p_ext_75v68  = p_ext_75v68_all(ext_idx);
p_ext_75vvpi = p_ext_75vvpi_all(ext_idx);
p_ext_68vvpi = p_ext_68vvpi_all(ext_idx);

%% 

function [p75v68, p75vvpi, p68vvpi] = run_ttests(d75, d68, dvpi)
    n = size(d75,2);
    p75v68  = zeros(1,n);
    p75vvpi = zeros(1,n);
    p68vvpi = zeros(1,n);
    for i = 1:n
        [~, p75v68(i)]  = ttest2(d75(:,i), d68(:,i),  'Vartype','unequal');
        [~, p75vvpi(i)] = ttest2(d75(:,i), dvpi(:,i), 'Vartype','unequal');
        [~, p68vvpi(i)] = ttest2(d68(:,i), dvpi(:,i), 'Vartype','unequal');
    end
end

p_int_75v68  = p_int_75v68_all(int_idx);
p_int_75vvpi = p_int_75vvpi_all(int_idx);
p_int_68vvpi = p_int_68vvpi_all(int_idx);

%% 

function s = star_label(p)
    if p < 0.001,    s = '***';
    elseif p < 0.01, s = '**';
    elseif p < 0.05, s = '*';
    else,            s = '';
    end
end

%% 
function plot_dotplot(d75, d68, dvpi, media_m, met_names, ...
        p75v68, p75vvpi, p68vvpi, alpha, col75, col68, colVPI, ttl)

    n   = length(met_names);
    jit = 0.15;

    Y75  = zeros(3,n); Y68 = zeros(3,n); YVPI = zeros(3,n);
    for i = 1:n
        Y75(:,i)  = d75(:,i)  - media_m(i);
        Y68(:,i)  = d68(:,i)  - media_m(i);
        YVPI(:,i) = dvpi(:,i) - media_m(i);
    end

    figure('Units','inches','Position',[1 1 max(8, n*1.2) 5.5]);
    ax = axes;
    hold(ax,'on');

    for i = 1:n
        y75  = Y75(:,i);
        y68  = Y68(:,i);
        yvpi = YVPI(:,i);

        xs75  = i - jit;
        xs68  = i;
        xsvpi = i + jit;

        scatter(ax, repmat(xs75, 3,1),  y75,  300, 'o', 'MarkerFaceColor',col75,  'MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
        scatter(ax, repmat(xs68, 3,1),  y68,  300, 'o', 'MarkerFaceColor',col68,  'MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
        scatter(ax, repmat(xsvpi,3,1),  yvpi, 300, 'o', 'MarkerFaceColor',colVPI, 'MarkerEdgeColor','none','MarkerFaceAlpha',0.85);

        lh = 0.08;
        plot(ax,[xs75-lh  xs75+lh],  [mean(y75)  mean(y75)],  '-','Color',col75,  'LineWidth',4.5);
        plot(ax,[xs68-lh  xs68+lh],  [mean(y68)  mean(y68)],  '-','Color',col68,  'LineWidth',4.5);
        plot(ax,[xsvpi-lh xsvpi+lh], [mean(yvpi) mean(yvpi)], '-','Color',colVPI, 'LineWidth',4.5);

        y_top = max([y75; y68; yvpi]);
        step  = (max(max([Y75; Y68; YVPI])) - min(min([Y75; Y68; YVPI]))) * 0.08;

        pairs  = {[xs75 xsvpi], [xs75 xs68], [xs68 xsvpi]};
        pvals  = [p75vvpi(i), p75v68(i), p68vvpi(i)];
        levels = [1, 2, 3];

        for k = 1:3
            s = star_label(pvals(k));
            if ~isempty(s)
                by = y_top + levels(k)*step;
                bx = pairs{k};
                plot(ax,[bx(1) bx(1)],[by-step*0.4 by],'k-','LineWidth',0.8);
                plot(ax,[bx(2) bx(2)],[by-step*0.4 by],'k-','LineWidth',0.8);
                plot(ax,[bx(1) bx(2)],[by by],          'k-','LineWidth',0.8);
                text(ax, mean(bx), by+step*0.15, s, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom',...
                    'FontSize',45,'FontWeight','bold','Color','k');
            end
        end
    end

    yline(ax, 0, '--k', 'LineWidth',0.8, 'Alpha',0.35);

    set(ax,'XTick',1:n,'XTickLabel',met_names,'XTickLabelRotation',35,...
        'FontSize',30,'Box','off','TickDir','out','XLim',[0.5 n+0.5]);
    ylabel(ax, 'log_2 fold change vs media','FontSize',20);
    title(ax, ttl,'FontSize',12,'FontWeight','bold');

    h(1) = scatter(ax,nan,nan,355,'o','MarkerFaceColor',col75, 'MarkerEdgeColor','none');
    h(2) = scatter(ax,nan,nan,355,'o','MarkerFaceColor',col68, 'MarkerEdgeColor','none');
    h(3) = scatter(ax,nan,nan,355,'o','MarkerFaceColor',colVPI,'MarkerEdgeColor','none');
    legend(ax,h,{'ST1-75','ST1-68','VPI10463'},'Location','northeast','FontSize',30,'Box','off');
end

%% 
compartments  = {'Extracellular', 'Intracellular'};
strain_data   = {st75_ext, st68_ext, vpi_ext; st75_int, st68_int, vpi_int};
media_means   = {media_ext_mean; media_int_mean};
met_names_all = {metabolites_ext; metabolites_int};

for c = 1:2
    fprintf('--- %s ---\n', compartments{c});
    fprintf('%-18s %-10s %-10s %-10s %-8s %-8s %-8s\n', ...
        'Metabolite', 'FC_75', 'FC_68', 'FC_VPI', 'p_75v68', 'p_75vVPI', 'p_68vVPI');
    fprintf('%s\n', repmat('-',1,75));

    d75  = strain_data{c,1};
    d68  = strain_data{c,2};
    dvpi = strain_data{c,3};
    mm   = media_means{c};
    mets = met_names_all{c};

    for i = 1:length(mets)
        
        fc75  = mean(d75(:,i))  - mm(i);
        fc68  = mean(d68(:,i))  - mm(i);
        fcvpi = mean(dvpi(:,i)) - mm(i);

        [~, p1] = ttest2(d75(:,i), d68(:,i),  'Vartype','unequal');
        [~, p2] = ttest2(d75(:,i), dvpi(:,i), 'Vartype','unequal');
        [~, p3] = ttest2(d68(:,i), dvpi(:,i), 'Vartype','unequal');

        sig = '';
        if p1 < 0.05 || p2 < 0.05 || p3 < 0.05
            sig = ' <-- SIG';
        end

        fprintf('%-18s %-10.2f %-10.2f %-10.2f %-8.3f %-8.3f %-8.3f%s\n', ...
            mets{i}, fc75, fc68, fcvpi, p1, p2, p3, sig);
    end
    fprintf('\n');
end

plot_dotplot(st75_ext(:,ext_idx), st68_ext(:,ext_idx), vpi_ext(:,ext_idx), ...
    media_ext_mean(ext_idx), metabolites_ext(ext_idx), ...
    p_ext_75v68, p_ext_75vvpi, p_ext_68vvpi, ...
    alpha, col75, col68, colVPI, ...
    'Extracellular metabolites (secretome)');

plot_dotplot(st75_int(:,int_idx), st68_int(:,int_idx), vpi_int(:,int_idx), ...
    media_int_mean(int_idx), metabolites_int(int_idx), ...
    p_int_75v68, p_int_75vvpi, p_int_68vvpi, ...
    alpha, col75, col68, colVPI, ...
    'Intracellular metabolites');
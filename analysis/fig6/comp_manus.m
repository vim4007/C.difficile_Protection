% data = readtable('/Users/vmishra/C.difficile_Protection/data/cerillo/competition.xlsx','NumHeaderLines', 0);
base_dir = fileparts(fileparts(mfilename('fullpath')));
data = readtable(fullfile(base_dir, 'data', 'cerillo', 'competition.xlsx'), 'NumHeaderLines', 0);

x    = double(string(data.Duration_Hours_));
function y = process_well(data, col,~)
    y = double(string(data.(col)));
    y = y - y(1);
    y(y < 0) = 0;
    y = smoothdata(y, 'loess', 40);
end

alanine_75_reps  = zeros(length(x), 3);
alanine_vpi_reps = zeros(length(x), 3);
for r = 1:3
    alanine_75_reps(:,r)  = process_well(data, sprintf('A%d', r), x);
    alanine_vpi_reps(:,r) = process_well(data, sprintf('B%d', r), x);
end

alanine_68_reps   = zeros(length(x), 3);
alanine_vpi2_reps = zeros(length(x), 3);
for r = 1:3
    alanine_68_reps(:,r)   = process_well(data, sprintf('A%d', r+3), x);
    alanine_vpi2_reps(:,r) = process_well(data, sprintf('B%d', r+3), x);
end

glycine_75_reps  = zeros(length(x), 3);
glycine_vpi_reps = zeros(length(x), 3);
for r = 1:3
    glycine_75_reps(:,r)  = process_well(data, sprintf('C%d', r), x);
    glycine_vpi_reps(:,r) = process_well(data, sprintf('D%d', r), x);
end

glycine_68_reps   = zeros(length(x), 3);
glycine_vpi2_reps = zeros(length(x), 3);
for r = 1:3
    glycine_68_reps(:,r)   = process_well(data, sprintf('C%d', r+3), x);
    glycine_vpi2_reps(:,r) = process_well(data, sprintf('D%d', r+3), x);
end

col75  = [0.20 0.53 0.74];
col68  = [0.20 0.63 0.29];
colVPI = [0.85 0.15 0.15];
%%
function plot_competition(ax, x, strain_reps, vpi_reps, strain_color, vpi_color, strain_name, ttl)
    hold(ax, 'on');

    sm  = mean(strain_reps, 2);
    ss  = std(strain_reps,  0, 2);
    vm  = mean(vpi_reps,    2);
    vs  = std(vpi_reps,     0, 2);

    t = x(:);

    % Shaded SD regions
    fill(ax, [t; flipud(t)], [sm+ss; flipud(sm-ss)], strain_color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill(ax, [t; flipud(t)], [vm+vs; flipud(vm-vs)], vpi_color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Mean lines
    plot(ax, t, sm, '-', 'Color', strain_color, 'LineWidth', 4.5, 'DisplayName', strain_name);
    plot(ax, t, vm, '-', 'Color', vpi_color,    'LineWidth', 4.5, 'DisplayName', 'VPI10463');

    set(ax, 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', ...
        'XLim', [0 max(t)], 'YLim', [0 1]);
    xlabel(ax, 'Time (hours)', 'FontSize', 20);
    ylabel(ax, 'OD_6_0_0_n_m', 'FontSize', 20);
    title(ax, ttl, 'FontSize', 15, 'FontWeight', 'bold');
    legend(ax, 'Location', 'northwest', 'FontSize', 22, 'Box', 'off');
end

%% 

figure('Units', 'inches', 'Position', [1 1 12 9]);

ax1 = subplot(2, 2, 1);
plot_competition(ax1, x, alanine_75_reps, alanine_vpi_reps, col75, colVPI, 'ST1-75', 'Alanine — ST1-75 vs VPI');

ax2 = subplot(2, 2, 2);
plot_competition(ax2, x, alanine_68_reps, alanine_vpi2_reps, col68, colVPI, 'ST1-68', 'Alanine — ST1-68 vs VPI');

ax3 = subplot(2, 2, 3);
plot_competition(ax3, x, glycine_75_reps, glycine_vpi_reps, col75, colVPI, 'ST1-75', 'Glycine — ST1-75 vs VPI');

ax4 = subplot(2, 2, 4);
plot_competition(ax4, x, glycine_68_reps, glycine_vpi2_reps, col68, colVPI, 'ST1-68', 'Glycine — ST1-68 vs VPI');

sgtitle('Competition Assay: Growth in Shared Medium', 'FontSize', 13, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');

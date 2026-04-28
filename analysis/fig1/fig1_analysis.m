%% Read the raw tables

% vir_tbl = readtable('/Users/vmishra/C.difficile_Protection/data/mouse/Scores/Virulence_screen_clean_table.csv');
% tbl = readtable('/Users/vmishra/C.difficile_Protection/data/mouse/Scores/ProtectionScreen_CDI_mouse.csv');

script_dir = fileparts(which('fig1_analysis'));
base_dir = fileparts(fileparts(script_dir));

vir_tbl = readtable(fullfile(base_dir, 'data', 'mouse', 'Scores', 'Virulence_screen_clean_table.csv'));
tbl = readtable(fullfile(base_dir, 'data', 'mouse', 'Scores', 'ProtectionScreen_CDI_mouse.csv'));
tbl = [tbl;vir_tbl];
tbl.cdiffstrain = string(tbl.cdiffstrain);
tbl.tx = string(tbl.tx);
tbl.experiment = string(tbl.experiment);
tbl.exp_id = tbl.experiment + "_" + tbl.cdiffstrain + "_" + tbl.mouse;
%% Make separate tables

s75secondary = tbl(endsWith(tbl.cdiffstrain, ".st1.75"), :); % co-infection with St1-75 as secondary infection
no_vpi_st175_idx = ~contains(tbl.cdiffstrain, ".vpi") & ~contains(tbl.cdiffstrain, ".st1.75");
primary = tbl(no_vpi_st175_idx, :); % "Virulence table from the protection screen = primary infection only"
vpi_ending_idx = endsWith(tbl.cdiffstrain, ".vpi");
vpi_alone_idx = strcmp(tbl.cdiffstrain, "vpi");
secondary = tbl(vpi_ending_idx | vpi_alone_idx, :);% "Protection table from the protection screen = co-infection with VPI"
%% Virulence 
primary.experiment = categorical(primary.experiment);
primary.exp_id = categorical(primary.exp_id);
primary.relweight(isnan(primary.relweight)) = 0;
primary.cdiffstrain = categorical(primary.cdiffstrain);
primary.cdiffstrain = reordercats(primary.cdiffstrain, ...
    ['ui'; setdiff(categories(primary.cdiffstrain), 'ui')]);
%%
primary.relweight = -1*primary.relweight;
virulence_model = fitlme(primary,'relweight ~ cdiffstrain + (1|day) + (1|exp_id)');

coeffs = fixedEffects(virulence_model);
names = virulence_model.CoefficientNames;
ci = coefCI(virulence_model);

% Combine into a table
virulence_score = table(names(:), coeffs(:), ci(:,1), ci(:,2), ...
    'VariableNames', {'cdiffstrain', 'Estimate', 'CI_Lower', 'CI_Upper'});
virulence_score = sortrows(virulence_score, 'Estimate', 'descend');
virulence_score.Strains = extractAfter(virulence_score.cdiffstrain, "cdiffstrain_");
virulence_score.Strains = string(virulence_score.Strains);

%% Protection
secondary.experiment = categorical(secondary.experiment);
secondary.exp_id = categorical(secondary.exp_id);
secondary.relweight(isnan(secondary.relweight)) = 0;
secondary.cdiffstrain = categorical(secondary.cdiffstrain);
secondary.cdiffstrain = reordercats(secondary.cdiffstrain, ...
    ['vpi'; setdiff(categories(secondary.cdiffstrain), 'vpi')]);
%%
protection_model = fitlme(secondary,'relweight ~ cdiffstrain + (1|day) + (1|exp_id)');

coeffs = fixedEffects(protection_model);
names = protection_model.CoefficientNames;
ci = coefCI(protection_model);

% Combine into a table
protection_score = table(names(:), coeffs(:), ci(:,1), ci(:,2), ...
    'VariableNames', {'cdiffstrain', 'Estimate', 'CI_Lower', 'CI_Upper'});
protection_score = sortrows(protection_score, 'Estimate', 'descend');
protection_score.Strains = extractAfter(protection_score.cdiffstrain, "cdiffstrain_");
protection_score.Strains = extractBefore(protection_score.Strains, ".vpi");
protection_score.Strains = string(protection_score.Strains);

%%
combined_score = protection_score.Strains;
combined_score = array2table(combined_score);
combined_score(1,:) = [];
combined_score = renamevars(combined_score, "combined_score", "Strains");

[lia, locb] = ismember(combined_score.Strains, virulence_score.Strains);
combined_score.Virulence_Estimate = NaN(height(combined_score), 1);
combined_score.Virulence_Estimate(lia) = virulence_score.Estimate(locb(lia));
combined_score.Virulence_CI_Lower = NaN(height(combined_score), 1);
combined_score.Virulence_CI_Upper = NaN(height(combined_score), 1);
combined_score.Virulence_CI_Lower(lia) = virulence_score.CI_Lower(locb(lia));
combined_score.Virulence_CI_Upper(lia) = virulence_score.CI_Upper(locb(lia));

[lia, locb] = ismember(combined_score.Strains, protection_score.Strains);
combined_score.Protection_Estimate = NaN(height(combined_score), 1);
combined_score.Protection_Estimate(lia) = protection_score.Estimate(locb(lia));
combined_score.Protection_CI_Lower = NaN(height(combined_score), 1);
combined_score.Protection_CI_Upper = NaN(height(combined_score), 1);
combined_score.Protection_CI_Lower(lia) = protection_score.CI_Lower(locb(lia));
combined_score.Protection_CI_Upper(lia) = protection_score.CI_Upper(locb(lia));

combined_score.Strains = upper(strrep(combined_score.Strains, '.', '-'));
combined_score(strcmp(combined_score.Strains, 'CD196'), :) = [];

%% Scatter with vir CIs

figure('Position', [100 100 800 600]);

hold on;
for i = 1:height(combined_score)
 
    virulence_err_lower = combined_score.Virulence_Estimate(i) - combined_score.Virulence_CI_Lower(i);
    virulence_err_upper = combined_score.Virulence_CI_Upper(i) - combined_score.Virulence_Estimate(i);
    
    protection_err_lower = combined_score.Protection_Estimate(i) - combined_score.Protection_CI_Lower(i);
    protection_err_upper = combined_score.Protection_CI_Upper(i) - combined_score.Protection_Estimate(i);
   
    errorbar(combined_score.Virulence_Estimate(i), combined_score.Protection_Estimate(i), ...
        protection_err_lower, protection_err_upper, virulence_err_lower, virulence_err_upper, ...
        'Color', [0.85 0.85 0.85], 'LineWidth', 1, 'CapSize', 3);
end

scatter(combined_score.Virulence_Estimate, combined_score.Protection_Estimate, ...
    100,'b', 'filled', 'MarkerFaceAlpha', 0.7);

mdl = fitlm(combined_score, 'Protection_Estimate ~ Virulence_Estimate');
x_range = linspace(min(combined_score.Virulence_Estimate), max(combined_score.Virulence_Estimate), 100);
newdata = table(x_range', 'VariableNames', {'Virulence_Estimate'});
y_fit = predict(mdl, newdata);
plot(x_range, y_fit, 'r-', 'LineWidth', 3);

text(combined_score.Virulence_Estimate - 1.3, combined_score.Protection_Estimate - 1.3, ...
    combined_score.Strains, 'FontSize', 22, 'Horiz', 'left', 'Vert', 'bottom');

xlabel('Virulence Score'); 
ylabel('Protection Score');
set(gca,'FontSize',20);

[R, P] = corrcoef(combined_score.Virulence_Estimate, combined_score.Protection_Estimate);
title(sprintf('(r = %.3f, p = %.5f, R^2 = %.3f)',R(1,2), P(1,2), R(1,2)^2),'FontSize',15);

hold off;
%%
primary.relweight = -1*primary.relweight;
%% Growth trajectories
strain_name = "st1.6.vpi";
secondary.cdiffstrain = string(secondary.cdiffstrain);
figure; 
hold on;
blue = [0 0.45 0.74];

strain_idx = strcmp(secondary.cdiffstrain, strain_name);
strain_data = secondary(strain_idx, :);

strain_data.exp_id = string(strain_data.exp_id);

all_exp_ids = unique(strain_data.exp_id);
died_exp_ids = [];
surviving_exp_ids = [];

for i = 1:numel(all_exp_ids)
    exp_id = all_exp_ids(i);
    mouse_data = strain_data(strain_data.exp_id == exp_id, :);

    if any(mouse_data.death == 1) || any(mouse_data.relweight == 0)
        died_exp_ids = [died_exp_ids; exp_id];
    else
        surviving_exp_ids = [surviving_exp_ids; exp_id];
    end
end

for i = 1:numel(died_exp_ids)
    exp_id = died_exp_ids(i);
    mouse_data = strain_data(strain_data.exp_id == exp_id, :);
    mouse_data = sortrows(mouse_data, 'day');

    non_zero_idx = find(mouse_data.relweight > 0, 1, 'last');
    
    if ~isempty(non_zero_idx)

        plot(mouse_data.day(1:non_zero_idx), mouse_data.relweight(1:non_zero_idx), ...
            'k-', 'LineWidth', 1);

        plot(mouse_data.day(non_zero_idx), mouse_data.relweight(non_zero_idx), ...
            'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

surviving_idx = ismember(strain_data.exp_id, surviving_exp_ids);
days_all = strain_data.day(surviving_idx);
exp_all = strain_data.exp_id(surviving_idx);
rel_all = strain_data.relweight(surviving_idx);
n_exps_surviving = numel(surviving_exp_ids);
n_exps_died = numel(died_exp_ids);

days = sort(unique(days_all));
mean_rel = NaN(numel(days), 1);
sd_rel = NaN(numel(days), 1);

for i = 1:numel(days)
    day_mask = days_all == days(i);
    day_exp = exp_all(day_mask);
    day_rel = rel_all(day_mask);

    if ~iscell(day_exp)
        day_exp = cellstr(string(day_exp));
    else
        day_exp = cellstr(day_exp);
    end
    
    unique_exps = unique(day_exp);
    exp_means = NaN(numel(unique_exps), 1);
    
    for j = 1:numel(unique_exps)
        exp_mask = strcmp(day_exp, unique_exps{j});
        exp_means(j) = nanmean(day_rel(exp_mask));
    end
    
    valid_means = exp_means(~isnan(exp_means));
    if ~isempty(valid_means)
        mean_rel(i) = nanmean(valid_means);
        if numel(valid_means) > 1
            sd_rel(i) = nanstd(valid_means);
        end
    end
end

valid = ~isnan(mean_rel);
x = days(valid); 
y = mean_rel(valid); 
e = sd_rel(valid);
e(isnan(e)) = 0; 

fill([x; flipud(x)], [y-e; flipud(y+e)], blue, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(x, y, '-', 'Color', blue, 'LineWidth', 4);

set(gca, 'FontSize', 40);
title(sprintf('%s (N=%d surviving, %d died)', strain_name, n_exps_surviving, n_exps_died), ...
    'FontSize', 40, 'FontWeight', 'bold');
xlabel('Day'); 

ylim([70 110]); 
xlim([0 7]);
grid on; 
box on;
hold off;
%% Bar plots for virulence and protection

combined_score_sorted = sortrows(combined_score, 'Virulence_Estimate', 'ascend');

err_lower = combined_score_sorted.Virulence_Estimate - combined_score_sorted.Virulence_CI_Lower;
err_upper = combined_score_sorted.Virulence_CI_Upper - combined_score_sorted.Virulence_Estimate;

figure;

bar(combined_score_sorted.Virulence_Estimate);
hold on;
errorbar(1:height(combined_score_sorted), combined_score_sorted.Virulence_Estimate, ...
    err_lower, err_upper, 'k.', 'LineWidth', 1.5);
hold off;
set(gca, 'XTick', 1:height(combined_score_sorted));
set(gca, 'XTickLabel', combined_score_sorted.Strains);
set(gca, 'XTickLabelRotation', 45);

xlabel('Strains', 'FontSize', 14);
ylabel('Virulence Estimate', 'FontSize', 14);ylim([-10 25]);
title('Virulence Estimates by Strain (with 95% CI)', 'FontSize', 16);
set(gca, 'FontSize', 20);
grid on;

set(gcf, 'Position', [100, 100, 1000, 600]);

combined_score_sorted = sortrows(combined_score, 'Protection_Estimate', 'ascend');

err_lower = combined_score_sorted.Protection_Estimate - combined_score_sorted.Protection_CI_Lower;
err_upper = combined_score_sorted.Protection_CI_Upper - combined_score_sorted.Protection_Estimate;

figure;
bar(combined_score_sorted.Protection_Estimate);
hold on;
errorbar(1:height(combined_score_sorted), combined_score_sorted.Protection_Estimate, ...
    err_lower, err_upper, 'k.', 'LineWidth', 1.5);
hold off;
set(gca, 'XTick', 1:height(combined_score_sorted));
set(gca, 'XTickLabel', combined_score_sorted.Strains);
set(gca, 'XTickLabelRotation', 45);
xlabel('Strains', 'FontSize', 14);
ylabel('Protection Estimate', 'FontSize', 14);ylim([-10 25]);
title('Protection Estimates by Strain (with 95% CI)', 'FontSize', 16);
set(gca, 'FontSize', 20);
grid on;
set(gcf, 'Position', [100, 100, 1000, 600]);
%% Survival curve

strains_of_interest = ["ui", "st1.75", "st1.2"];
filtered_data = primary(ismember(primary.cdiffstrain, strains_of_interest), :);

figure;
hold on;
colors = {'b', 'r', 'g'};
linestyles = {'-', '--', ':'}; 
for i = 1:length(strains_of_interest)
    strain = strains_of_interest(i);
    strain_data = filtered_data(filtered_data.cdiffstrain == strain, :);
    if isempty(strain_data)
        warning('No data found for strain: %s', strain);
        continue;
    end
    [G, exp_ids] = findgroups(strain_data.exp_id);
    final_time = splitapply(@max, strain_data.day, G);
    final_death = splitapply(@max, strain_data.death, G);
    
    if sum(final_death) == 0
        x = [0; sort(unique(final_time)); max(final_time)*1.1];
        f = ones(size(x));
    else
        [f, x] = ecdf(final_time, 'Censoring', ~final_death, 'Function', 'survivor');
        
        max_time = max(final_time);
        if x(end) < max_time
            x = [x; max_time];
            f = [f; f(end)];
        end
    end

    stairs(x, f, 'Color', colors{i}, 'LineStyle', linestyles{i}, ...
        'LineWidth', 5.5, 'DisplayName', ...
        sprintf('%s', strain));
end

xlabel('Time (days)', 'FontSize', 12);
ylabel('Survival Probability', 'FontSize', 12);

set(gca,'FontSize',25);
legend('Location', 'southwest', 'FontSize', 20);
grid on;
ylim([0 1.05]);
xlim([0 inf]);
hold off;

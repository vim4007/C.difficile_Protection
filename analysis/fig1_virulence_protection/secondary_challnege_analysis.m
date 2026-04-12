tbl = readtable('weights.xlsx');
%%
% Step 1: Create Day column (0-indexed from first date)
% Convert to datetime if not already
tbl.Date = datetime(tbl.Date);

% Calculate days from first date
first_date = min(tbl.Date);
tbl.Day = days(tbl.Date - first_date);

% Step 2: Create Group column from ExperimentalGroup
tbl.Group = cell(height(tbl), 1);
for i = 1:height(tbl)
    first_letter = tbl.ExperimentalGroup{i}(1);
    switch first_letter
        case 'A'
            tbl.Group{i} = 'ST1-75';
        case 'B'
            tbl.Group{i} = 'ST1-68';
        case 'C'
            tbl.Group{i} = 'VPI';
    end
end

% Step 3: Create relweight column (normalized to day 0 = 100 for each mouse)
tbl.relweight = nan(height(tbl), 1);
unique_mice = unique(tbl.ExperimentalGroup);

for i = 1:length(unique_mice)
    mouse_idx = strcmp(tbl.ExperimentalGroup, unique_mice{i});
    mouse_data = tbl(mouse_idx, :);
    
    % Find day 0 weight for this mouse
    day0_weight = mouse_data.Weight_g_(mouse_data.Day == 0);
    
    % Normalize all weights for this mouse
    tbl.relweight(mouse_idx) = (mouse_data.Weight_g_ / day0_weight) * 100;
end

% Create death column (1 where relweight is NaN, 0 otherwise)
tbl.death = double(isnan(tbl.relweight));

% Replace NaN relweight values with 0
tbl.relweight(isnan(tbl.relweight)) = 0;

tbl.death = (tbl.relweight == 0);
tbl.death = double(tbl.death);
%% Plot ST1-75
% Filter for ST1-75 group only
st1_75_data = tbl(strcmp(tbl.Group, 'ST1-75'), :);

% Get unique days
unique_days = unique(st1_75_data.Day);
mean_relweight = nan(length(unique_days), 1);
std_relweight = nan(length(unique_days), 1);

% Calculate mean and std for each day
for i = 1:length(unique_days)
    day_data = st1_75_data(st1_75_data.Day == unique_days(i), :);
    mean_relweight(i) = mean(day_data.relweight);
    std_relweight(i) = std(day_data.relweight);
end

% Plot
figure;
hold on;

% Shaded area (mean ± std)
upper = mean_relweight + std_relweight;
lower = mean_relweight - std_relweight;
fill([unique_days; flipud(unique_days)], [upper; flipud(lower)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Mean line
plot(unique_days, mean_relweight, 'b-', 'LineWidth', 2);
xline(6, 'k--', 'LineWidth', 1.5);
xline(23, 'k--', 'LineWidth', 1.5);
xlabel('Day', 'FontSize', 14);ylim([80 115]);xlim([0 34]);
ylabel('Relative Weight', 'FontSize', 14);
set(gca,'FontSize',30);
title('ST1-75', 'FontSize', 30);
grid on;
hold off;


%% Plot ST1-68
st1_68_data = tbl(strcmp(tbl.Group, 'ST1-68'), :);
unique_mice_68 = unique(st1_68_data.ExperimentalGroup);
unique_days_68 = sort(unique(st1_68_data.Day));

% Separate surviving and dead mice
died_68 = {};
surviving_68 = {};
for m = 1:length(unique_mice_68)
    mouse_data = st1_68_data(strcmp(st1_68_data.ExperimentalGroup, unique_mice_68{m}), :);
    mouse_data = sortrows(mouse_data, 'Day');
    if any(mouse_data.death == 1)
        died_68{end+1} = unique_mice_68{m};
    else
        surviving_68{end+1} = unique_mice_68{m};
    end
end

% Calculate mean and SD from surviving mice only
mean_relweight_68 = nan(length(unique_days_68), 1);
std_relweight_68 = nan(length(unique_days_68), 1);
for i = 1:length(unique_days_68)
    day_vals = [];
    for m = 1:length(surviving_68)
        mouse_data = st1_68_data(strcmp(st1_68_data.ExperimentalGroup, surviving_68{m}), :);
        day_idx = mouse_data.Day == unique_days_68(i);
        if any(day_idx)
            day_vals = [day_vals; mouse_data.relweight(day_idx)];
        end
    end
    if ~isempty(day_vals)
        mean_relweight_68(i) = mean(day_vals);
        if length(day_vals) > 1
            std_relweight_68(i) = std(day_vals);
        else
            std_relweight_68(i) = 0;
        end
    end
end

% Plot
figure;
hold on;
blue = [0 0.45 0.74];

% Plot individual dead mouse traces
for m = 1:length(died_68)
    mouse_data = st1_68_data(strcmp(st1_68_data.ExperimentalGroup, died_68{m}), :);
    mouse_data = sortrows(mouse_data, 'Day');
    death_idx = find(mouse_data.death == 1, 1, 'first');
    if death_idx > 1
        days_plot = mouse_data.Day(1:death_idx);
        weights_plot = mouse_data.relweight(1:death_idx-1);
        weights_plot = [weights_plot; weights_plot(end)];
        plot(days_plot, weights_plot, 'k-', 'LineWidth', 1);
        plot(mouse_data.Day(death_idx), weights_plot(end), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Shaded SD and mean line from survivors
valid = ~isnan(mean_relweight_68);
x = unique_days_68(valid);
y = mean_relweight_68(valid);
e = std_relweight_68(valid);
e(isnan(e)) = 0;
fill([x; flipud(x)], [y-e; flipud(y+e)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x, y, 'b-', 'LineWidth', 2);

xline(6, 'k--', 'LineWidth', 1.5);
xline(23, 'k--', 'LineWidth', 1.5);
xlabel('Day', 'FontSize', 14);
ylabel('Relative Weight', 'FontSize', 14);
ylim([80 115]); xlim([0 34]);
set(gca, 'FontSize', 30);
title('ST1-68', 'FontSize', 30);
grid on;
hold off;

%% Plot VPI
vpi_data = tbl(strcmp(tbl.Group, 'VPI'), :);
unique_mice_vpi = unique(vpi_data.ExperimentalGroup);
unique_days_vpi = sort(unique(vpi_data.Day));

% Separate surviving and dead mice
died_vpi = {};
surviving_vpi = {};
for m = 1:length(unique_mice_vpi)
    mouse_data = vpi_data(strcmp(vpi_data.ExperimentalGroup, unique_mice_vpi{m}), :);
    mouse_data = sortrows(mouse_data, 'Day');
    if any(mouse_data.death == 1)
        died_vpi{end+1} = unique_mice_vpi{m};
    else
        surviving_vpi{end+1} = unique_mice_vpi{m};
    end
end

% Calculate mean and SD from all mice up to death day
mean_relweight_vpi = nan(length(unique_days_vpi), 1);
std_relweight_vpi = nan(length(unique_days_vpi), 1);
for i = 1:length(unique_days_vpi)
    day = unique_days_vpi(i);
    day_vals = [];
    for m = 1:length(unique_mice_vpi)
        mouse_data = vpi_data(strcmp(vpi_data.ExperimentalGroup, unique_mice_vpi{m}), :);
        mouse_data = sortrows(mouse_data, 'Day');
        death_idx = find(mouse_data.death == 1, 1, 'first');
        if isempty(death_idx)
            % Never died - include if day exists
            day_idx = mouse_data.Day == day;
            if any(day_idx)
                day_vals = [day_vals; mouse_data.relweight(day_idx)];
            end
        else
            % Only include if this day is before death
            death_day = mouse_data.Day(death_idx);
            if day < death_day
                day_idx = mouse_data.Day == day;
                if any(day_idx)
                    day_vals = [day_vals; mouse_data.relweight(day_idx)];
                end
            end
        end
    end
    if ~isempty(day_vals)
        mean_relweight_vpi(i) = mean(day_vals);
        if length(day_vals) > 1
            std_relweight_vpi(i) = std(day_vals);
        else
            std_relweight_vpi(i) = 0;
        end
    end
end

% Plot
figure;
hold on;
blue = [0 0.45 0.74];

% Plot individual dead mouse traces
for m = 1:length(died_vpi)
    mouse_data = vpi_data(strcmp(vpi_data.ExperimentalGroup, died_vpi{m}), :);
    mouse_data = sortrows(mouse_data, 'Day');
    death_idx = find(mouse_data.death == 1, 1, 'first');
    if death_idx > 1
        days_plot = mouse_data.Day(1:death_idx);
        weights_plot = mouse_data.relweight(1:death_idx-1);
        weights_plot = [weights_plot; weights_plot(end)];
        plot(days_plot, weights_plot, 'k-', 'LineWidth', 1);
        plot(mouse_data.Day(death_idx), weights_plot(end), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Shaded SD and mean line from survivors
valid = ~isnan(mean_relweight_vpi);
x = unique_days_vpi(valid);
y = mean_relweight_vpi(valid);
e = std_relweight_vpi(valid);
e(isnan(e)) = 0;
fill([x; flipud(x)], [y-e; flipud(y+e)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x, y, 'b-', 'LineWidth', 2);

xline(6, 'k--', 'LineWidth', 1.5);
xline(23, 'k--', 'LineWidth', 1.5);
xlabel('Day', 'FontSize', 14);
ylabel('Relative Weight', 'FontSize', 14);
ylim([80 115]); xlim([0 34]);
set(gca, 'FontSize', 30);
title('VPI10463', 'FontSize', 30);
grid on;
hold off;

%% Survival

% Survival curves for all 3 groups
figure;
hold on;

groups = {'ST1-75', 'ST1-68', 'VPI'};
colors = {'r', 'g', 'b'};

for g = 1:length(groups)
    group_data = tbl(strcmp(tbl.Group, groups{g}), :);
    
    % Get unique mice and days
    unique_mice = unique(group_data.ExperimentalGroup);
    unique_days = sort(unique(group_data.Day));
    
    % Calculate survival proportion at each day
    survival_prop = nan(length(unique_days), 1);
    
    for i = 1:length(unique_days)
        day = unique_days(i);
        n_alive = 0;
        
        for m = 1:length(unique_mice)
            mouse_data = group_data(strcmp(group_data.ExperimentalGroup, unique_mice{m}), :);
            mouse_data = sortrows(mouse_data, 'Day');
            
            % Find first death day
            death_idx = find(mouse_data.death == 1, 1, 'first');
            
            if isempty(death_idx)
                % Never died - still alive
                n_alive = n_alive + 1;
            else
                death_day = mouse_data.Day(death_idx);
                if day < death_day
                    % Still alive on this day
                    n_alive = n_alive + 1;
                end
            end
        end
        
        survival_prop(i) = n_alive / length(unique_mice) * 100;
    end
    
    % Plot survival curve for this group
    stairs(unique_days, survival_prop, colors{g}, 'LineWidth', 2, 'DisplayName', groups{g});
end

xlabel('Day', 'FontSize', 14);
ylabel('Survival (%)', 'FontSize', 14);
title('Survival Curves by Group', 'FontSize', 15);
xlim([0 35]);
ylim([0 105]);
legend('Location', 'southwest');
grid on;
set(gca, 'FontSize', 20);
hold off;
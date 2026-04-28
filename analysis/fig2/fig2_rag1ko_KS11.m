

%tbl = readtable('/Users/vmishra/C.difficile_Protection/data/mouse/rag1ko/KS11_rechallenge_relweight.xlsx');
script_dir = fileparts(which('fig2_rag1ko_KS11'));
base_dir = fileparts(fileparts(script_dir));
tbl = readtable(fullfile(base_dir, 'data', 'mouse', 'rag1ko', 'KS11_rechallenge_relweight.xlsx'));
tbl.group = string(tbl.group);
tbl.mouse_id = string(tbl.mouse_id);

groupKeys  = {'uninfected_b6',   'st175_b6',    'st175_rag1ko'};
groupTitles = {'Uninfected B6', 'ST1.75 B6', 'ST1.75 RAG1 KO'};
blue = [0 0.45 0.74];

%% 
for g = 1:3
    grp = tbl(tbl.group == groupKeys{g}, :);
    mouse_ids = unique(grp.mouse_id);
    nMice = length(mouse_ids);

    died_mice = string([]);
    surviving_mice = string([]);

    for m = 1:nMice
        mid = mouse_ids(m);
        mdata = grp(grp.mouse_id == mid, :);
        postDay0 = mdata.relweight(mdata.day > 0);
        if any(postDay0 == 0 | isnan(postDay0))
            died_mice = [died_mice; mid];
        else
            surviving_mice = [surviving_mice; mid];
        end
    end

    figure;
    hold on;
    for m = 1:length(died_mice)
        mid = died_mice(m);
        mdata = sortrows(grp(grp.mouse_id == mid, :), 'day');
        traj = mdata.relweight;
        days = mdata.day;
        lastAlive = find(traj > 0 & ~isnan(traj), 1, 'last');
        if ~isempty(lastAlive)
            plot(days(1:lastAlive), traj(1:lastAlive), ...
                'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
            plot(days(lastAlive), traj(lastAlive), ...
                'rx', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    if ~isempty(surviving_mice)
        all_days = sort(unique(grp.day));
        mean_rel = NaN(length(all_days), 1);
        sd_rel   = NaN(length(all_days), 1);

        for d = 1:length(all_days)
            day_vals = [];
            for m = 1:length(surviving_mice)
                mid = surviving_mice(m);
                mdata = grp(grp.mouse_id == mid & grp.day == all_days(d), :);
                if ~isempty(mdata) && ~isnan(mdata.relweight)
                    day_vals = [day_vals; mdata.relweight];
                end
            end
            if ~isempty(day_vals)
                mean_rel(d) = mean(day_vals, 'omitnan');
                if length(day_vals) > 1
                    sd_rel(d) = std(day_vals, 'omitnan');
                else
                    sd_rel(d) = 0;
                end
            end
        end

        sd_rel(isnan(sd_rel)) = 0;
        x = all_days;
        y = mean_rel;
        e = sd_rel;

        fill([x; flipud(x)], [y-e; flipud(y+e)], blue, ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(x, y, '-', 'Color', blue, 'LineWidth', 4);
    end

    xlabel('Day post secondary challenge');
    ylabel('Relative weight (%)');
    title(groupTitles{g});
    ylim([70 115]); xlim([0 7]);
    set(gca, 'FontSize', 20);
    grid on; box on;

    fprintf('Group: %s | Surviving: %d | Died: %d\n', ...
        groupTitles{g}, length(surviving_mice), length(died_mice));

    hold off;
end

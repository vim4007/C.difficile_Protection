weights = readtable('/Users/vmishra/C.difficile_Protection/data/qPCR/weights.xlsx');
qpcr = readtable("/Users/vmishra/C.difficile_Protection/data/qPCR/qPCR section.xlsx");
weights = renamevars(weights, "Var1", "Days");
weights(1,:) = [];
%% --- Weight trajectories ---
days_w = weights.Days;

wt_st175 = weights.ST1_75;
wt_vpi   = weights.VPI;
wt_mix   = [weights.mix, weights.mix_1, weights.mix_2];  % 3x3

mix_mean = mean(wt_mix, 2);
mix_std  = std(wt_mix, 0, 2);

figure('Position',[100 100 600 450]);
hold on;
fill([days_w; flipud(days_w)], ...
     [mix_mean - mix_std; flipud(mix_mean + mix_std)], ...
     [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
     'HandleVisibility','off');


plot(days_w, mix_mean, '-o', 'Color',[0.4 0.4 0.4], ...
    'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',[0.4 0.4 0.4], ...
    'DisplayName','ST1-75 + VPI');

plot(days_w, wt_st175, '-o', 'Color',[0.20 0.55 0.80], ...
    'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',[0.20 0.55 0.80], ...
    'DisplayName','ST1-75');

plot(days_w, wt_vpi, '-o', 'Color',[0.85 0.33 0.10], ...
    'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',[0.85 0.33 0.10], ...
    'DisplayName','VPI');

yline(100,'k--','LineWidth',1,'HandleVisibility','off');
xlabel('Days post-infection');
ylabel('Relative Weight');
title('Weight trajectories');
legend('Location','southwest','Box','off','FontSize',18);
xlim([0 3.2]); ylim([60 120]);
set(gca,'XTick',[0 2 3],'FontSize',20);
box off;

%% --- qPCR stacked bar plots ---

days_q   = [0.5, 1, 3];
col_blue = [0.0  0.0  0.8];
col_red  = [1  0.0  0.0];

days_q = [0.5, 1, 3];

for d = 1:3
    idx = qpcr.Day == days_q(d);
    day_data = qpcr(idx,:);
    
    st175_alone(d) = day_data.ST1_75(day_data.Mouse == 1);

    vpi_alone(d) = day_data.VPI(day_data.Mouse == 2);

    co_idx = ismember(day_data.Mouse, [3,4,5]);
    mean_co_st175(d) = nanmean(day_data.ST1_75(co_idx));
    mean_co_vpi(d)   = nanmean(day_data.VPI(co_idx));
end

figure('Position',[100 100 1200 450]);

subplot(1,3,1);
for d = 1:3
    bar(d, st175_alone(d), 'FaceColor', col_blue, 'EdgeColor','k','LineWidth',1.2);
    hold on;
end
set(gca,'XTick',1:3,'XTickLabel',{'0.5','1','3'},'FontSize',14,'FontWeight','bold');
xlabel('Day','FontSize',16,'FontWeight','bold');
ylabel('ST1-75 (%)','FontSize',16,'FontWeight','bold');
ylim([0 100]); xlim([0.4 3.6]);
set(gca,'YGrid','on','GridColor',[0.8 0.8 0.8],'GridAlpha',1);
box on; set(gca,'LineWidth',1.2);

subplot(1,3,2);
for d = 1:3
    if ~isnan(vpi_alone(d))
        bar(d, vpi_alone(d), 'FaceColor', col_red, 'EdgeColor','k','LineWidth',1.2);
    else
        bar(d, 100, 'FaceColor','none','EdgeColor','k','LineWidth',1.2);
    end
    hold on;
end
set(gca,'XTick',1:3,'XTickLabel',{'0.5','1','3'},'FontSize',14,'FontWeight','bold');
xlabel('Day','FontSize',16,'FontWeight','bold');
ylabel('VPI (%)','FontSize',16,'FontWeight','bold');
ylim([0 100]); xlim([0.4 3.6]);
set(gca,'YGrid','on','GridColor',[0.8 0.8 0.8],'GridAlpha',1);
box on; set(gca,'LineWidth',1.2);

subplot(1,3,3);
b = bar(1:3, [mean_co_vpi', mean_co_st175'], 'stacked', ...
    'EdgeColor','k','LineWidth',1.2);
b(1).FaceColor = col_red;
b(2).FaceColor = col_blue;

set(gca,'XTick',1:3,'XTickLabel',{'0.5','1','3'},'FontSize',14,'FontWeight','bold');
xlabel('Day','FontSize',16,'FontWeight','bold');
ylabel('Percentage','FontSize',16,'FontWeight','bold');
ylim([0 100]); xlim([0.4 3.6]);
set(gca,'YGrid','on','GridColor',[0.8 0.8 0.8],'GridAlpha',1);
legend({'VPI','ST1-75'},'Location','northwest','Box','on','FontSize',12);
box on; set(gca,'LineWidth',1.2);
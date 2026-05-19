%% Parameters

p_common.e     = 0.5;      
p_common.d     = 0.05;     
p_common.D     = 0.2;     
p_common.delta0 = 0.08;
p_common.us     = 0.25;
p_common.u      = 0.25;     
p_common.U     = 0.06;     

tspan = [0 24];            

Rs0      = 10.0;            
R1_0     = 10.0;            
R2_0     = 10.0;          

%% ----- Pair 1: ST1-75 (strain 1) vs VPI10463 (strain 2) ---------------
pairA = p_common;
pairA.n1 = 50;  pairA.n2 = 6;  pairA.ns = 27;
pairA.label  = 'ST1-75 vs VPI10463';
pairA.name1  = 'ST1-75';
pairA.name2  = 'VPI10463';

%% ----- Pair 2: ST1-68 (strain 1) vs VPI10463 (strain 2) ---------------
pairB = p_common;
pairB.U  = 0;              
pairB.n1 = 11; pairB.n2 = 21; pairB.ns = 12;
pairB.label  = 'ST1-68 vs VPI10463';
pairB.name1  = 'ST1-68';
pairB.name2  = 'VPI10463';

%% ----- Integrate ------------------------------------------------------
y0 = [1; 1; Rs0; R1_0; R2_0];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

[tA, yA] = ode45(@(t,y) cr_rhs(t,y,pairA), tspan, y0, opts);
[tB, yB] = ode45(@(t,y) cr_rhs(t,y,pairB), tspan, y0, opts);

%% 
fig1 = figure('Units','inches','Position',[1 1 12 7],'Color','w');

col_str1 = [0.20 0.45 0.85];   
col_str2 = [0.85 0.15 0.15];   
col_Rs   = [0.00 0.00 0.00];   
col_R1   = [0.20 0.45 0.85];   
col_R2   = [0.85 0.15 0.15];   

subplot(2,2,1); hold on;
plot(tA, yA(:,1),'-','Color',col_str1,'LineWidth',4.2);
plot(tA, yA(:,2),'-','Color',col_str2,'LineWidth',4.2);
ylabel('Population density'); xlabel('Time (hours)');ylim([0 15]);
title(pairA.label,'FontWeight','bold');
legend({pairA.name1, pairA.name2},'Location','best','Box','off');
grid on; box off; set(gca,'FontSize',20);

subplot(2,2,3); hold on;
plot(tA, yA(:,3),'-', 'Color',col_Rs,'LineWidth',4.0);
plot(tA, yA(:,4),'--','Color',col_R1,'LineWidth',4.0);
plot(tA, yA(:,5),'--','Color',col_R2,'LineWidth',4.0);
ylabel('Resource concentration'); xlabel('Time (hours)');ylim([0 5]);
legend({'R_s (shared)','R_1 (private 1)','R_2 (private 2)'}, ...
       'Location','best','Box','off');
grid on; box off; set(gca,'FontSize',20);

subplot(2,2,2); hold on;
plot(tB, yB(:,1),'-','Color',[0.20 0.63 0.29],'LineWidth',4.2);  % green for ST1-68
plot(tB, yB(:,2),'-','Color',col_str2,'LineWidth',4.2);
ylabel('Population density'); xlabel('Time (hours)');ylim([0 15]);
title(pairB.label,'FontWeight','bold');
legend({pairB.name1, pairB.name2},'Location','best','Box','off');
grid on; box off; set(gca,'FontSize',20);

subplot(2,2,4); hold on;
plot(tB, yB(:,3),'-', 'Color',col_Rs,'LineWidth',4.0);
plot(tB, yB(:,4),'--','Color',[0.20 0.63 0.29],'LineWidth',4.0);
plot(tB, yB(:,5),'--','Color',col_R2,'LineWidth',4.0);
ylabel('Resource concentration'); xlabel('Time (hours)');ylim([0 5]);
legend({'R_s (shared)','R_1 (private 1)','R_2 (private 2)'}, ...
       'Location','best','Box','off');
grid on; box off; set(gca,'FontSize',20);

%%
function dy = cr_rhs(~, y, p)

N1 = y(1); N2 = y(2); Rs = y(3); R1 = y(4); R2 = y(5);

uf = p.us*(1 + p.U);  
dtilde = p.d + p.D;

dN1 = N1*(p.e*uf*Rs + p.e*p.u*R1 - dtilde);
dN2 = N2*(p.e*p.us*Rs + p.e*p.u*R2 - dtilde);

dRs = p.ns*p.delta0 - p.D*Rs - Rs*(uf*N1 + p.us*N2);
dR1 = p.n1*p.delta0 - p.D*R1 - p.u*N1*R1;
dR2 = p.n2*p.delta0 - p.D*R2 - p.u*N2*R2;

dy = [dN1; dN2; dRs; dR1; dR2];
end


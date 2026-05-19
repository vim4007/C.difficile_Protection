%% CR Model Phase Portraits 

e      = 0.5;
d      = 0.05;
D      = 0.2;
delta0 = 0.08;
u      = 0.25;
dtilde = d + D;
tspan = [0 150];
Rs0 = 10; R1_0 = 10; R2_0 = 10;
opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'NonNegative',1:5);
fs_small = 20;
fs_large = 25;
%% Pair definitions
pairs(1).n1    = 50; pairs(1).n2 = 6;  pairs(1).ns = 27;
pairs(1).name1 = 'ST1-75';
pairs(1).name2 = 'VPI10463';
pairs(1).label = 'ST1-75 vs VPI10463';
pairs(1).col1  = [0.20 0.45 0.85];   % blue
pairs(1).xlim  = 14;
pairs(1).ylim  = 13;
pairs(1).fname = 'cr_portrait_ST175';

pairs(2).n1    = 11; pairs(2).n2 = 21; pairs(2).ns = 12;
pairs(2).name1 = 'ST1-68';
pairs(2).name2 = 'VPI10463';
pairs(2).label = 'ST1-68 vs VPI10463';
pairs(2).col1  = [0.20 0.63 0.29];   % green
pairs(2).xlim  = 7;
pairs(2).ylim  = 13;
pairs(2).fname = 'cr_portrait_ST168';

col_VPI = [0.85 0.15 0.15];
col_eq  = [0.25 0.25 0.25];
col_5x  = [0.65 0.10 0.80];


cmap4 = [0.86 0.91 1.00;   % blue  -- strain1 only
         1.00 0.88 0.88;   % red   -- VPI only
         0.88 1.00 0.88;   % green -- both grow
         0.94 0.94 0.94];  % gray  -- neither

%% Loop over pairs -- one figure each
for pk = 1:2

    p = pairs(pk);
    xl = p.xlim; yl = p.ylim;

    fig = figure('Units','inches','Position',[1 1 10 8],'Color','w');
    set(fig,'Renderer','painters');
    ax = axes('Parent',fig); hold on;

    %% 1. Shaded growth regions (original contourf approach)
    Ng  = 250;
    N1g = linspace(0, xl, Ng);
    N2g = linspace(0, yl, Ng);
    G1  = zeros(Ng);
    G2  = zeros(Ng);
    for i = 1:Ng
        for j = 1:Ng
            N1v  = N1g(j); N2v = N2g(i);
            Rs_v = p.ns*delta0 / (D + u*(N1v + N2v + 1e-9));
            R1_v = p.n1*delta0 / (D + u*(N1v + 1e-9));
            R2_v = p.n2*delta0 / (D + u*(N2v + 1e-9));
            G1(i,j) = e*u*Rs_v + e*u*R1_v - dtilde;
            G2(i,j) = e*u*Rs_v + e*u*R2_v - dtilde;
        end
    end
    Z = zeros(Ng);
    Z((G1 > 0) & (G2 < 0)) = 1;   % strain1 only
    Z((G1 < 0) & (G2 > 0)) = 2;   % VPI only
    Z((G1 > 0) & (G2 > 0)) = 3;   % both grow
    Z((G1 < 0) & (G2 < 0)) = 4;   % neither
    contourf(N1g, N2g, Z, [0.5 1.5 2.5 3.5 4.5], 'LineStyle','none');
    colormap(ax, cmap4); caxis([0.5 4.5]);

    %% 2. Vector field
    Nq  = 13;
    N1q = linspace(0.01, xl*0.98, Nq);
    N2q = linspace(0.01, yl*0.98, Nq);
    [QN1, QN2] = meshgrid(N1q, N2q);
    QdN1 = zeros(Nq); QdN2 = zeros(Nq);
    for i = 1:Nq
        for j = 1:Nq
            N1v       = QN1(i,j); N2v = QN2(i,j);
            Rs_v      = p.ns*delta0 / (D + u*(N1v + N2v));
            R1_v      = p.n1*delta0 / (D + u*N1v);
            R2_v      = p.n2*delta0 / (D + u*N2v);
            QdN1(i,j) = N1v*(e*u*Rs_v + e*u*R1_v - dtilde);
            QdN2(i,j) = N2v*(e*u*Rs_v + e*u*R2_v - dtilde);
        end
    end
    mag = sqrt(QdN1.^2 + QdN2.^2) + 1e-12;
    quiver(QN1, QN2, QdN1./mag, QdN2./mag, 0.45, ...
        'Color',[0.55 0.55 0.55],'LineWidth',0.9,'MaxHeadSize',0.5);

    %% 3. Nullclines
    N2v = linspace(0.001, yl*1.02, 800);
    NC1 = nan(size(N2v));
    NC2 = nan(size(N2v));
    for i = 1:length(N2v)
        n2 = N2v(i);
        % N1 nullcline: G1(N1, n2) = 0
        f1 = @(N1) e*u*(p.ns*delta0/(D + u*(N1 + n2))) ...
                 + e*u*(p.n1*delta0/(D + u*N1)) - dtilde;
        try
            s = fzero(f1, [1e-4, xl*4]);
            if s > 0 && s <= xl*1.05; NC1(i) = s; end
        catch; end
        % N2 nullcline: G2(N1, n2) = 0, solve for N1
        R2_v = p.n2*delta0 / (D + u*n2);
        rhs  = dtilde - e*u*R2_v;
        if rhs > 0
            N1_sol = (e*u*p.ns*delta0/rhs - D)/u - n2;
            if N1_sol > 0 && N1_sol <= xl*1.05; NC2(i) = N1_sol; end
        end
    end
    hnc1 = plot(NC1, N2v, '-', 'Color',p.col1,  'LineWidth',3.5);
    hnc2 = plot(NC2, N2v, '-', 'Color',col_VPI, 'LineWidth',3.5);

    %% 4. Trajectories
    y0_eq = [1; 1; Rs0; R1_0; R2_0];
    y0_5x = [1; 5; Rs0; R1_0; R2_0];
    [~, y_eq] = ode45(@(t,y) cr_rhs(t,y,p,e,u,d,D,delta0), tspan, y0_eq, opts);
    [~, y_5x] = ode45(@(t,y) cr_rhs(t,y,p,e,u,d,D,delta0), tspan, y0_5x, opts);

    hteq = plot(y_eq(:,1), y_eq(:,2), '-',  'Color',col_eq, 'LineWidth',2.8);
    ht5x = plot(y_5x(:,1), y_5x(:,2), '--', 'Color',col_5x, 'LineWidth',3.2);

    draw_arrow(y_eq, 0.35, col_eq, ax);
    draw_arrow(y_5x, 0.25, col_5x, ax);

    % Starting points
    scatter(1, 1, 150, col_eq,'filled','o','MarkerEdgeColor','k','LineWidth',1.2);
    scatter(1, 5, 190, col_5x,'filled','p','MarkerEdgeColor','k','LineWidth',1.5);
    text(1.2, 5.4, '\bf VPI_0 = 5\times','FontSize',fs_small, ...
        'Color',col_5x,'FontWeight','bold');

    %% 5. Steady states

    % E0
    scatter(0, 0, 130, [0.5 0.5 0.5],'filled','o','MarkerEdgeColor','k');
    text(0.08*xl, -0.06*yl, 'E_0','FontSize',fs_small, ...
        'Color',[0.4 0.4 0.4],'FontWeight','bold');

    % E1: strain1 alone (N2=0)
    f_E1 = @(N1) e*u*(p.ns*delta0/(D + u*N1)) ...
                + e*u*(p.n1*delta0/(D + u*N1)) - dtilde;
    N1E1 = fzero(f_E1, [0.1, xl*4]);
    scatter(N1E1, 0, 230, p.col1,'o','MarkerEdgeColor',p.col1,'LineWidth',2.5);
    text(N1E1, 0.05*yl, sprintf('E_1\n(%.1f, 0)',N1E1), ...
        'FontSize',fs_small,'Color',p.col1,'FontWeight','bold', ...
        'HorizontalAlignment','center');
    text(N1E1, -0.09*yl, '[unstable]','FontSize',fs_small-4,'Color',p.col1, ...
        'HorizontalAlignment','center');

    % E2: VPI alone (N1=0)
    f_E2 = @(N2) e*u*(p.ns*delta0/(D + u*N2)) ...
                + e*u*(p.n2*delta0/(D + u*N2)) - dtilde;
    N2E2 = fzero(f_E2, [0.1, yl*4]);
    scatter(0, N2E2, 230, col_VPI,'o','MarkerEdgeColor',col_VPI,'LineWidth',2.5);
    text(0.06*xl, N2E2+0.05*yl, sprintf('E_2\n(0, %.1f)',N2E2), ...
        'FontSize',fs_small,'Color',col_VPI,'FontWeight','bold');
    text(0.06*xl, N2E2-0.09*yl, '[unstable]','FontSize',fs_small-4,'Color',col_VPI);

    % E*: coexistence
    Estar = find_coexistence(p, e, u, d, D, delta0, dtilde, xl, yl);
    if ~isempty(Estar)
        scatter(Estar(1), Estar(2), 290, [0.5 0.2 0.7],'filled','d', ...
            'MarkerEdgeColor','k','LineWidth',1.5);
        text(Estar(1)+0.04*xl, Estar(2)+0.04*yl, ...
            sprintf('E^* [stable]\n(%.2f, %.2f)',Estar(1),Estar(2)), ...
            'FontSize',fs_small,'Color',[0.5 0.2 0.7],'FontWeight','bold');
    end

    %% 6. Parameter annotation (top-left)
    mu1 = e*u*(p.ns*delta0/D) + e*u*(p.n1*delta0/D) - dtilde;
    mu2 = e*u*(p.ns*delta0/D) + e*u*(p.n2*delta0/D) - dtilde;
    ann = sprintf('n_1=%d, n_2=%d, n_s=%d\n\\mu_1 = %+.3f\n\\mu_2 = %+.3f', ...
        p.n1, p.n2, p.ns, mu1, mu2);
    text(xl*0.02, yl*0.98, ann,'FontSize',fs_small,'VerticalAlignment','top', ...
        'BackgroundColor','w','EdgeColor',[0.5 0.5 0.5],'Margin',5);

    %% 7. Axes
    xlabel(sprintf('N_1  (%s)',p.name1),'FontSize',fs_small);
    ylabel('N_2  (VPI10463)',           'FontSize',fs_small);
    title(p.label,'FontWeight','bold',  'FontSize',fs_large);
    xlim([0 xl]); ylim([0 yl]);
    set(ax,'FontSize',fs_small,'Layer','top'); box on; grid off;

    %% 8. Legend -- placed outside to the right, no overlap
    leg = legend([hnc1, hnc2, hteq, ht5x], ...
        {sprintf('N_1 nullcline (%s)',p.name1), ...
         'N_2 nullcline (VPI10463)', ...
         'Trajectory: N_1=N_2=1', ...
         'Trajectory: N_2=5 (VPI 5\times)'}, ...
        'Box','on','FontSize',fs_small);
    leg.Location = 'none';
    leg.Units    = 'normalized';
    % Shrink axes slightly to make room, then place legend to the right
    ax.Position(3) = 0.58;   % reduce axes width
    drawnow;
    axpos          = ax.Position;
    leg.Position(1) = axpos(1) + axpos(3) + 0.02;
    leg.Position(2) = axpos(2) + (axpos(4) - leg.Position(4))*0.5;

    

end  % end pair loop

%% ---- Helper functions ------------------------------------------------------

function dy = cr_rhs(~, y, p, e, u, d, D, delta0)
    N1 = max(y(1),0); N2 = max(y(2),0);
    Rs = max(y(3),0); R1 = max(y(4),0); R2 = max(y(5),0);
    dtilde = d + D;
    dN1 = N1*(e*u*Rs + e*u*R1 - dtilde);
    dN2 = N2*(e*u*Rs + e*u*R2 - dtilde);
    dRs = p.ns*delta0 - D*Rs - u*Rs*(N1 + N2);
    dR1 = p.n1*delta0 - D*R1 - u*N1*R1;
    dR2 = p.n2*delta0 - D*R2 - u*N2*R2;
    dy  = [dN1; dN2; dRs; dR1; dR2];
end

function Estar = find_coexistence(p, e, u, d, D, delta0, dtilde, xl, yl)
    F = @(x) [ e*u*(p.ns*delta0/(D + u*(x(1)+x(2)))) ...
                + e*u*(p.n1*delta0/(D + u*x(1))) - dtilde; ...
               e*u*(p.ns*delta0/(D + u*(x(1)+x(2)))) ...
                + e*u*(p.n2*delta0/(D + u*x(2))) - dtilde ];
    fsopts = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);
    try
        [sol, ~, flag] = fsolve(F, [xl*0.3; yl*0.3], fsopts);
        if flag > 0 && sol(1) > 0.01 && sol(2) > 0.01 ...
                    && sol(1) < xl*2   && sol(2) < yl*2
            Estar = sol;
        else
            Estar = [];
        end
    catch
        Estar = [];
    end
end

function draw_arrow(y_traj, frac, col, ax)
    n   = size(y_traj,1);
    idx = max(1, round(n*frac));
    if idx >= n; idx = n-1; end
    x  = y_traj(idx,  1);   y  = y_traj(idx,  2);
    dx = y_traj(idx+1,1)-x; dy = y_traj(idx+1,2)-y;
    len = sqrt(dx^2+dy^2)+1e-12;
    dx  = dx/len; dy = dy/len;
    sc  = 0.015*(ax.XLim(2)-ax.XLim(1));
    xd  = [x-sc*(dx+0.4*dy),  x,  x-sc*(dx-0.4*dy)];
    yd  = [y-sc*(dy-0.4*dx),  y,  y-sc*(dy+0.4*dx)];
    plot(ax, xd, yd, '-','Color',col,'LineWidth',2.5);
end
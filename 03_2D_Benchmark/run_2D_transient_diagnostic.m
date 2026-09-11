clear; clc; close all;

Lx_m = 200;
Ly_m = 100;

Lx = 201;
Ly = 101;

dx = Lx_m/(Lx-1);
dy = Ly_m/(Ly-1);
dt = 0.1;
e = dx/dt;

tau = 0.55;

u0 = 0.01;
v0 = 0.00;
D = 0.02;

C0 = 1.0;
x0 = 10.0;
y0 = 50.0;

ix0 = round(x0/dx) + 1;
iy0 = round(y0/dy) + 1;

times = [10; 100; 500; 1000; 2000; 4000];
iters = round(times/dt);

ns = length(times);
Nt = round(max(times)/dt);

neg_tol = -1e-12;

solid = zeros(Lx,Ly);
Cen = zeros(Lx,Ly);

u = u0*ones(Lx,Ly);
v = v0*ones(Lx,Ly);

Cen(ix0,iy0) = C0;

[ex,ey] = setup(e);

feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
f = feq;

min_feq0 = min(feq(:));
neg_feq0 = sum(feq(:) < neg_tol);

xc = (0:Lx-1)'*dx;
yc = (0:Ly-1)*dy;
[X,Y] = ndgrid(xc,yc);

minC = zeros(ns,1);
maxC = zeros(ns,1);
negN = zeros(ns,1);
massV = zeros(ns,1);

cx = zeros(ns,1);
cy = zeros(ns,1);

peakX = zeros(ns,1);
peakY = zeros(ns,1);
peakC = zeros(ns,1);

sigX = zeros(ns,1);
sigY = zeros(ns,1);

expX = zeros(ns,1);
expY = zeros(ns,1);
sigTheory = zeros(ns,1);

hasNeg = false(ns,1);

fields = zeros(Lx,Ly,ns);

hist_step = 100;
nh = floor(Nt/hist_step);

hist_t = zeros(nh,1);
hist_min = zeros(nh,1);
hist_neg = zeros(nh,1);
hist_mass = zeros(nh,1);

hi = 0;

first_neg = NaN;
first_recover = NaN;
neg_seen = false;

fprintf('2D transient diagnostic\n');
fprintf('Domain = %.0f x %.0f m\n', Lx_m,Ly_m);
fprintf('Grid = %d x %d\n', Lx,Ly);
fprintf('dx = %.4f m, dy = %.4f m, dt = %.4f s\n', dx,dy,dt);
fprintf('u = %.4f m/s, v = %.4f m/s\n', u0,v0);
fprintf('D = %.4f m^2/s, tau = %.2f\n', D,tau);
fprintf('Source = (%.2f, %.2f) m\n', x0,y0);
fprintf('Initial min feq = %.12e\n', min_feq0);
fprintf('Initial negative populations = %d\n', neg_feq0);
fprintf('Iterations = %d\n\n', Nt);

si = 1;

for iter = 1:Nt

    ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);
    ftemp = Peri_BC_inOut(1,Ly,ftemp,Lx);
    ftemp = Peri_BC_norSou(1,Lx,ftemp,Ly);

    [Cen,f] = solution(ftemp,1,Lx,1,Ly);

    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

    tc = iter*dt;
    cur_min = min(Cen(:));
    cur_neg = sum(Cen(:) < neg_tol);

    if cur_neg > 0
        neg_seen = true;

        if isnan(first_neg)
            first_neg = tc;
            fprintf('First negative C at %.2f s, min C = %.12e\n', ...
                first_neg,cur_min);
        end
    end

    if neg_seen && isnan(first_recover) && cur_neg == 0
        first_recover = tc;
    end

    if mod(iter,hist_step) == 0
        hi = hi + 1;

        hist_t(hi) = tc;
        hist_min(hi) = cur_min;
        hist_neg(hi) = cur_neg;
        hist_mass(hi) = sum(Cen(:))*dx*dy;
    end

    if si <= ns && iter == iters(si)

        tc = times(si);

        cmin = min(Cen(:));
        cmax = max(Cen(:));
        nneg = sum(Cen(:) < neg_tol);
        mass = sum(Cen(:))*dx*dy;

        [pC,ind] = max(Cen(:));
        [pxi,pyi] = ind2sub(size(Cen),ind);

        px = (pxi-1)*dx;
        py = (pyi-1)*dy;

        c_x = sum(X(:).*Cen(:))*dx*dy/mass;
        c_y = sum(Y(:).*Cen(:))*dx*dy/mass;

        sx = sqrt(sum(((X(:)-c_x).^2).*Cen(:))*dx*dy/mass);
        sy = sqrt(sum(((Y(:)-c_y).^2).*Cen(:))*dx*dy/mass);

        expt_x = x0 + u0*tc;
        expt_y = y0 + v0*tc;
        st = sqrt(2*D*tc);

        minC(si) = cmin;
        maxC(si) = cmax;
        negN(si) = nneg;
        massV(si) = mass;

        peakX(si) = px;
        peakY(si) = py;
        peakC(si) = pC;

        cx(si) = c_x;
        cy(si) = c_y;

        sigX(si) = sx;
        sigY(si) = sy;

        expX(si) = expt_x;
        expY(si) = expt_y;
        sigTheory(si) = st;

        hasNeg(si) = nneg > 0;
        fields(:,:,si) = Cen;

        fprintf('\nt = %.0f s\n', tc);
        fprintf('Min C = %.12e, Max C = %.10f\n', cmin,cmax);
        fprintf('Negative points = %d\n', nneg);
        fprintf('Mass = %.10f\n', mass);
        fprintf('Peak = (%.2f, %.2f) m\n', px,py);
        fprintf('Centroid = (%.6f, %.6f) m\n', c_x,c_y);
        fprintf('Expected = (%.6f, %.6f) m\n', expt_x,expt_y);
        fprintf('Sigma x = %.6f, Sigma y = %.6f, theory = %.6f\n', ...
            sx,sy,st);

        si = si + 1;
    end
end

Time_s = times;

summary_table = table( ...
    Time_s,minC,maxC,negN,massV,peakX,peakY,cx,expX,cy,expY, ...
    sigX,sigY,sigTheory,hasNeg, ...
    'VariableNames',{ ...
    'Time_s','MinimumC','MaximumC','NegativePointCount','TotalMass', ...
    'PeakX','PeakY','CentroidX','ExpectedX','CentroidY','ExpectedY', ...
    'SigmaX','SigmaY','TheoreticalSigma','HasNegativeConcentration'});

fprintf('\n');
disp(summary_table);

if isnan(first_neg)
    fprintf('No meaningful negative concentration detected.\n');
else
    fprintf('First negative concentration = %.2f s\n', first_neg);
end

if isnan(first_recover)
    fprintf('Negative concentrations did not disappear by %.0f s.\n', max(times));
else
    fprintf('First later time with no negative points = %.2f s\n', ...
        first_recover);
end

project_folder = fileparts(fileparts(mfilename('fullpath')));
resdir = fullfile(project_folder,'06_Results','2D_Transient_Diagnostic');

if ~exist(resdir,'dir')
    mkdir(resdir);
end

writetable(summary_table, ...
    fullfile(resdir,'2D_transient_diagnostic_summary.csv'));

save(fullfile(resdir,'2D_transient_diagnostic.mat'));

figure;
plot(hist_t(1:hi),hist_min(1:hi),'LineWidth',1.5);
hold on;
yline(neg_tol,'--','Negative threshold');
xlabel('Time (s)');
ylabel('Minimum concentration');
title('2D Early-Time Minimum Concentration');
grid on;
hold off;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_transient_minimum_concentration.png'), ...
    'Resolution',300);

figure;
plot(hist_t(1:hi),hist_neg(1:hi),'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Number of negative grid points');
title('2D Negative Concentration Points Through Time');
grid on;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_transient_negative_count.png'), ...
    'Resolution',300);

figure;
tiledlayout(2,3);

for k = 1:ns

    nexttile;

    contourf(xc,yc,fields(:,:,k)',20);
    colorbar;

    centre_x = x0 + u0*times(k);

    xlim([max(0,centre_x-45) min(Lx_m,centre_x+45)]);
    ylim([10 90]);

    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('t = %.0f s',times(k)));
end

sgtitle('2D Point-Source Transient Evolution');

exportgraphics(gcf, ...
    fullfile(resdir,'2D_transient_snapshots.png'), ...
    'Resolution',300);

fprintf('\nDiagnostic completed.\n');
fprintf('Results saved in: %s\n',resdir);
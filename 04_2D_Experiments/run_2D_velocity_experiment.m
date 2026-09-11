clear; clc; close all;

exp_dir = fileparts(mfilename('fullpath'));
proj_dir = fileparts(exp_dir);

bm_dir = fullfile(proj_dir,'03_2D_Benchmark');
old_dir = fullfile(proj_dir,'01_1D_Benchmark');

if contains(path,old_dir)
    rmpath(old_dir);
end

addpath(bm_dir,'-begin');
clear functions;
rehash path;

Lxm = 200;
Lym = 100;

Lx = 201;
Ly = 101;

dx = 1;
dy = 1;
dt = 0.1;
e = dx/dt;

tau = 0.55;
D = 0.020;
v0 = 0.0;

uvals = [0.0050; 0.0075; 0.0100; 0.0125; 0.0150];
nc = length(uvals);

C0 = 1.0;
x0 = 10.0;
y0 = 50.0;

ix0 = round(x0/dx) + 1;
iy0 = round(y0/dy) + 1;

tfinal = 4000;
Nt = round(tfinal/dt);

xu = (0:Lx-2)'*dx;
yu = (0:Ly-2)*dy;

[X,Y] = ndgrid(xu,yu);

Nx = length(xu);
Ny = length(yu);

Pe_v = zeros(nc,1);

rmse_v = zeros(nc,1);
mae_v = zeros(nc,1);
r2_v = zeros(nc,1);
maxerr_v = zeros(nc,1);

cmin_v = zeros(nc,1);
neg_v = zeros(nc,1);

mass_num = zeros(nc,1);
mass_ana = zeros(nc,1);

pk_num = zeros(nc,1);
pk_ana = zeros(nc,1);

pkx_num = zeros(nc,1);
pky_num = zeros(nc,1);
pkx_ana = zeros(nc,1);
pky_ana = zeros(nc,1);

cx_v = zeros(nc,1);
cy_v = zeros(nc,1);

xexp_v = zeros(nc,1);
yexp_v = zeros(nc,1);
cxerr_v = zeros(nc,1);

sx_v = zeros(nc,1);
sy_v = zeros(nc,1);
sig_ref = zeros(nc,1);

hasNeg = false(nc,1);
runtime_v = zeros(nc,1);

Cnum_all = zeros(Nx,Ny,nc);
Cana_all = zeros(Nx,Ny,nc);

fprintf('2D velocity experiment: %d cases, t = %.0f s\n\n',nc,tfinal);

for c = 1:nc

    timer = tic;

    u0 = uvals(c);
    Pe = u0*dx/D;

    Pe_v(c) = Pe;

    xexp = mod(x0 + u0*tfinal,Lxm);
    yexp = mod(y0 + v0*tfinal,Lym);

    xexp_v(c) = xexp;
    yexp_v(c) = yexp;

    solid = zeros(Lx,Ly);
    Cen = zeros(Lx,Ly);

    u = u0*ones(Lx,Ly);
    v = v0*ones(Lx,Ly);

    Cen(ix0,iy0) = C0;

    [ex,ey] = setup(e);

    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
    f = feq;

    fprintf('u = %.4f, Pe = %.4f, min feq = %.4e\n', ...
        u0,Pe,min(feq(:)));

    for iter = 1:Nt

        ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);
        ftemp = Peri_BC_inOut(1,Ly,ftemp,Lx);
        ftemp = Peri_BC_norSou(1,Lx,ftemp,Ly);

        [Cen,f] = solution(ftemp,1,Lx,1,Ly);

        feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
    end

    Cnum = Cen(1:Lx-1,1:Ly-1);

    Cana = zeros(size(X));
    A = C0/(4*pi*D*tfinal);

    for mx = -1:1
        for my = -1:1

            xi = xexp + mx*Lxm;
            yi = yexp + my*Lym;

            Cana = Cana + A .* exp( ...
                -((X-xi).^2 + (Y-yi).^2)/(4*D*tfinal));
        end
    end

    err = Cnum - Cana;

    rmse = sqrt(mean(err(:).^2));
    mae = mean(abs(err(:)));

    SS_res = sum(err(:).^2);
    SS_tot = sum((Cana(:)-mean(Cana(:))).^2);

    r2 = 1 - SS_res/SS_tot;
    maxerr = max(abs(err(:)));

    cmin = min(Cnum(:));
    nneg = sum(Cnum(:) < -1e-12);

    mnum = sum(Cnum(:))*dx*dy;
    mana = sum(Cana(:))*dx*dy;

    [pnum,idx] = max(Cnum(:));
    [ix,iy] = ind2sub(size(Cnum),idx);

    pnx = xu(ix);
    pny = yu(iy);

    [pana,idx] = max(Cana(:));
    [ix,iy] = ind2sub(size(Cana),idx);

    pax = xu(ix);
    pay = yu(iy);

    thx = 2*pi*X/Lxm;

    ax = atan2( ...
        sum(Cnum(:).*sin(thx(:))), ...
        sum(Cnum(:).*cos(thx(:))) );

    if ax < 0
        ax = ax + 2*pi;
    end

    cx = Lxm*ax/(2*pi);

    thy = 2*pi*Y/Lym;

    ay = atan2( ...
        sum(Cnum(:).*sin(thy(:))), ...
        sum(Cnum(:).*cos(thy(:))) );

    if ay < 0
        ay = ay + 2*pi;
    end

    cy = Lym*ay/(2*pi);

    dcx = abs(cx-xexp);
    cxerr = min(dcx,Lxm-dcx);

    rx = mod(X-cx+Lxm/2,Lxm)-Lxm/2;
    ry = mod(Y-cy+Lym/2,Lym)-Lym/2;

    sx = sqrt(sum((rx(:).^2).*Cnum(:))*dx*dy/mnum);
    sy = sqrt(sum((ry(:).^2).*Cnum(:))*dx*dy/mnum);

    st = sqrt(2*D*tfinal);

    rmse_v(c) = rmse;
    mae_v(c) = mae;
    r2_v(c) = r2;
    maxerr_v(c) = maxerr;

    cmin_v(c) = cmin;
    neg_v(c) = nneg;

    mass_num(c) = mnum;
    mass_ana(c) = mana;

    pk_num(c) = pnum;
    pk_ana(c) = pana;

    pkx_num(c) = pnx;
    pky_num(c) = pny;

    pkx_ana(c) = pax;
    pky_ana(c) = pay;

    cx_v(c) = cx;
    cy_v(c) = cy;

    cxerr_v(c) = cxerr;

    sx_v(c) = sx;
    sy_v(c) = sy;
    sig_ref(c) = st;

    hasNeg(c) = nneg > 0;

    Cnum_all(:,:,c) = Cnum;
    Cana_all(:,:,c) = Cana;

    runtime_v(c) = toc(timer);

    fprintf('  RMSE = %.4e, R2 = %.8f, negatives = %d\n', ...
        rmse,r2,nneg);

    fprintf('  mass: num = %.6f, ana = %.6f\n',mnum,mana);

    fprintf('  peak = %.6f at (%.1f, %.1f), centre = (%.3f, %.3f)\n', ...
        pnum,pnx,pny,cx,cy);

    fprintf('  expected x = %.3f, x error = %.4f\n',xexp,cxerr);

    fprintf('  sigma x = %.4f, sigma y = %.4f, theory = %.4f\n\n', ...
        sx,sy,st);
end

T = table( ...
    uvals,Pe_v,rmse_v,mae_v,r2_v,maxerr_v,cmin_v,neg_v, ...
    mass_num,mass_ana,pk_num,pk_ana,pkx_num,pky_num, ...
    pkx_ana,pky_ana,cx_v,cy_v,xexp_v,yexp_v,cxerr_v, ...
    sx_v,sy_v,sig_ref,hasNeg,runtime_v, ...
    'VariableNames',{ ...
    'Velocity','Peclet','RMSE','MAE','R2','MaxAbsoluteError', ...
    'MinimumC','NegativePointCount','NumericalMass','AnalyticalMass', ...
    'NumericalPeak','AnalyticalPeak','NumericalPeakX','NumericalPeakY', ...
    'AnalyticalPeakX','AnalyticalPeakY','CentroidX','CentroidY', ...
    'ExpectedCentreX','ExpectedCentreY','CentroidErrorX', ...
    'SigmaX','SigmaY','TheoreticalSigma','HasNegativeConcentration', ...
    'Runtime_seconds'});

fprintf('\n');
disp(T);

resdir = fullfile(proj_dir,'06_Results','2D_Velocity');

if ~exist(resdir,'dir')
    mkdir(resdir);
end

writetable(T,fullfile(resdir,'2D_velocity_summary.csv'));
save(fullfile(resdir,'2D_velocity_results.mat'));

figure;
tiledlayout(2,3);

for c = 1:nc

    nexttile;

    contourf(xu,yu,Cnum_all(:,:,c)',20);
    hold on;

    contour(xu,yu,Cana_all(:,:,c)',8,'k--','LineWidth',1);

    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');

    title(sprintf('u = %.4f m/s, Pe = %.3f',uvals(c),Pe_v(c)));

    axis equal;
    xlim([0 100]);
    ylim([0 100]);

    hold off;
end

sgtitle('Effect of Velocity on 2D Pollutant Transport at t = 4000 s');

exportgraphics(gcf, ...
    fullfile(resdir,'2D_velocity_concentration_fields.png'), ...
    'Resolution',300);

figure;

plot(uvals,cx_v,'o-','LineWidth',1.5);
hold on;

plot(uvals,xexp_v,'s--','LineWidth',1.5);

xlabel('Velocity u (m/s)');
ylabel('Plume centre x-position (m)');
title('2D Plume Position vs Flow Velocity');

legend('LBM Numerical','Expected x_0 + ut','Location','best');

grid on;
hold off;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_velocity_centroid_vs_velocity.png'), ...
    'Resolution',300);

figure;

plot(uvals,sx_v,'o-','LineWidth',1.5);
hold on;

plot(uvals,sy_v,'s-','LineWidth',1.5);
plot(uvals,sig_ref,'--','LineWidth',1.5);

xlabel('Velocity u (m/s)');
ylabel('Plume standard deviation (m)');
title('2D Plume Width vs Flow Velocity');

legend( ...
    'Numerical \sigma_x', ...
    'Numerical \sigma_y', ...
    'Theoretical \sigma', ...
    'Location','best');

grid on;
hold off;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_velocity_sigma_vs_velocity.png'), ...
    'Resolution',300);

figure;

plot(uvals,pk_num,'o-','LineWidth',1.5);
hold on;

plot(uvals,pk_ana,'s--','LineWidth',1.5);

xlabel('Velocity u (m/s)');
ylabel('Peak concentration');
title('Peak Concentration vs Flow Velocity');

legend('LBM Numerical','Analytical','Location','best');

grid on;
hold off;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_velocity_peak_vs_velocity.png'), ...
    'Resolution',300);

figure;

plot(uvals,r2_v,'o-','LineWidth',1.5);

xlabel('Velocity u (m/s)');
ylabel('R^2');
title('2D Analytical Validation Across Velocity Values');

grid on;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_velocity_R2.png'), ...
    'Resolution',300);

fprintf('Done, results in %s\n',resdir);
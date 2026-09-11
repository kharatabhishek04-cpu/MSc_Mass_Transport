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

Lch = 200;
Lx = 201;

dx = 1;
dy = 1;
dt = 0.1;
e = dx/dt;

tau = 0.55;

u0 = 0.01;
v0 = 0.0;
D = 0.020;

widths = [20; 30; 40; 60; 100];
nc = length(widths);

C0 = 1.0;
x0 = 10.0;
ix0 = round(x0/dx) + 1;

tfinal = 4000;
Nt = round(tfinal/dt);

xexp = mod(x0 + u0*tfinal,Lch);
sig0 = sqrt(2*D*tfinal);

gridLy = zeros(nc,1);
srcY = zeros(nc,1);

cmin_v = zeros(nc,1);
neg_v = zeros(nc,1);
mass_v = zeros(nc,1);

pk_v = zeros(nc,1);
pkx_v = zeros(nc,1);
pky_v = zeros(nc,1);

cx_v = zeros(nc,1);
cy_v = zeros(nc,1);
cxerr_v = zeros(nc,1);

sx_v = zeros(nc,1);
sy_v = zeros(nc,1);
syratio_v = zeros(nc,1);

mean_v = zeros(nc,1);
runtime_v = zeros(nc,1);
hasNeg_v = false(nc,1);

fields = cell(nc,1);
ys = cell(nc,1);

fprintf('Channel width experiment: %d cases, t = %.0f s\n\n',nc,tfinal);

for c = 1:nc

    timer = tic;

    W = widths(c);
    Ly = round(W/dy) + 1;

    gridLy(c) = Ly;

    y0 = W/2;
    srcY(c) = y0;
    iy0 = round(y0/dy) + 1;

    xu = (0:Lx-2)'*dx;
    yc = (0:Ly-1)*dy;

    [X,Y] = ndgrid(xu,yc);

    solid = zeros(Lx,Ly);
    Cen = zeros(Lx,Ly);

    u = u0*ones(Lx,Ly);
    v = v0*ones(Lx,Ly);

    Cen(ix0,iy0) = C0;

    [ex,ey] = setup(e);

    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
    f = feq;

    fprintf('Width %.0f m, grid %d x %d, min feq %.4e\n', ...
        W,Lx,Ly,min(feq(:)));

    for iter = 1:Nt

        ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);
        ftemp = Peri_BC_inOut(1,Ly,ftemp,Lx);

        for x = 1:Lx

            fn = f(2,x,Ly) - (f(2,x,Ly)-feq(2,x,Ly))/tau;
            ftemp(4,x,Ly) = fn;

            fs = f(4,x,1) - (f(4,x,1)-feq(4,x,1))/tau;
            ftemp(2,x,1) = fs;
        end

        [Cen,f] = solution(ftemp,1,Lx,1,Ly);

        feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
    end

    C = Cen(1:Lx-1,:);

    cmin = min(C(:));
    nneg = sum(C(:) < -1e-12);
    mass = sum(C(:))*dx*dy;

    [pk,idx] = max(C(:));
    [ix,iy] = ind2sub(size(C),idx);

    pkx = xu(ix);
    pky = yc(iy);

    thx = 2*pi*X/Lch;

    ax = atan2( ...
        sum(C(:).*sin(thx(:))), ...
        sum(C(:).*cos(thx(:))) );

    if ax < 0
        ax = ax + 2*pi;
    end

    cx = Lch*ax/(2*pi);
    cy = sum(Y(:).*C(:))*dx*dy/mass;

    dxc = abs(cx-xexp);
    cxerr = min(dxc,Lch-dxc);

    rx = mod(X-cx+Lch/2,Lch)-Lch/2;
    ry = Y-cy;

    sx = sqrt(sum((rx(:).^2).*C(:))*dx*dy/mass);
    sy = sqrt(sum((ry(:).^2).*C(:))*dx*dy/mass);

    syratio = sy/sig0;
    cmean = mean(C(:));

    cmin_v(c) = cmin;
    neg_v(c) = nneg;
    mass_v(c) = mass;

    pk_v(c) = pk;
    pkx_v(c) = pkx;
    pky_v(c) = pky;

    cx_v(c) = cx;
    cy_v(c) = cy;
    cxerr_v(c) = cxerr;

    sx_v(c) = sx;
    sy_v(c) = sy;
    syratio_v(c) = syratio;

    mean_v(c) = cmean;
    hasNeg_v(c) = nneg > 0;

    fields{c} = C;
    ys{c} = yc;

    runtime_v(c) = toc(timer);

    fprintf('  min C = %.4e, negatives = %d, mass = %.6f\n', ...
        cmin,nneg,mass);

    fprintf('  peak = %.6f at (%.1f, %.1f), centre = (%.3f, %.3f)\n', ...
        pk,pkx,pky,cx,cy);

    fprintf('  sigma x = %.4f, sigma y = %.4f, ratio = %.4f\n\n', ...
        sx,sy,syratio);
end

expX = xexp*ones(nc,1);
sigRef = sig0*ones(nc,1);

T = table( ...
    widths,gridLy,srcY,cmin_v,neg_v,mass_v,pk_v,pkx_v,pky_v, ...
    cx_v,cy_v,expX,cxerr_v,sx_v,sy_v,sigRef,syratio_v, ...
    mean_v,hasNeg_v,runtime_v, ...
    'VariableNames',{ ...
    'ChannelWidth','GridLy','SourceY','MinimumC','NegativePointCount', ...
    'TotalMass','NumericalPeak','PeakX','PeakY','CentroidX','CentroidY', ...
    'ExpectedCentreX','CentroidErrorX','SigmaX','SigmaY','UnboundedSigma', ...
    'SigmaYRatio','MeanConcentration','HasNegativeConcentration', ...
    'Runtime_seconds'});

fprintf('\n');
disp(T);

resdir = fullfile(proj_dir,'06_Results','2D_Channel_Width');

if ~exist(resdir,'dir')
    mkdir(resdir);
end

writetable(T,fullfile(resdir,'2D_channel_width_summary.csv'));
save(fullfile(resdir,'2D_channel_width_results.mat'));

figure;
tiledlayout(2,3);

for c = 1:nc

    nexttile;

    contourf(xu,ys{c},fields{c}',20);
    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('Width = %.0f m',widths(c)));

    xlim([0 100]);
    ylim([0 widths(c)]);

    axis equal;
end

sgtitle('Effect of Channel Width on 2D Pollutant Transport at t = 4000 s');

exportgraphics(gcf, ...
    fullfile(resdir,'2D_channel_width_concentration_fields.png'), ...
    'Resolution',300);

figure;

plot(widths,pk_v,'o-','LineWidth',1.5);

xlabel('Channel width (m)');
ylabel('Peak concentration');
title('Peak Concentration vs Channel Width');

grid on;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_channel_width_peak.png'), ...
    'Resolution',300);

figure;

plot(widths,sy_v,'o-','LineWidth',1.5);
hold on;

yline(sig0,'--','Unbounded \sigma');

xlabel('Channel width (m)');
ylabel('\sigma_y (m)');
title('Lateral Plume Spread vs Channel Width');

grid on;
hold off;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_channel_width_sigma_y.png'), ...
    'Resolution',300);

figure;

plot(widths,sx_v,'o-','LineWidth',1.5);
hold on;

yline(sig0,'--','Unbounded \sigma');

xlabel('Channel width (m)');
ylabel('\sigma_x (m)');
title('Longitudinal Plume Spread vs Channel Width');

grid on;
hold off;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_channel_width_sigma_x.png'), ...
    'Resolution',300);

figure;

plot(widths,cx_v,'o-','LineWidth',1.5);
hold on;

yline(xexp,'--','Expected x centre');

xlabel('Channel width (m)');
ylabel('Plume centroid x (m)');
title('Plume Transport Position vs Channel Width');

grid on;
hold off;

exportgraphics(gcf, ...
    fullfile(resdir,'2D_channel_width_centroid_x.png'), ...
    'Resolution',300);

fprintf('Done, results in %s\n',resdir);
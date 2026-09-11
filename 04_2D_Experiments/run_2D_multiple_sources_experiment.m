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

u0 = 0.01;
v0 = 0.0;
D = 0.020;

sx0 = [10; 25; 40];
sy0 = [25; 50; 75];
smass = [1; 1; 1];

nsrc = length(sx0);
mass0 = sum(smass);

tsnap = [2000; 3000; 4000];
iter_snap = round(tsnap/dt);

Ns = length(tsnap);
Nt = round(max(tsnap)/dt);

xu = (0:Lx-2)'*dx;
yu = (0:Ly-2)*dy;

[X,Y] = ndgrid(xu,yu);

Nx = length(xu);
Ny = length(yu);

solid = zeros(Lx,Ly);
Cen = zeros(Lx,Ly);

u = u0*ones(Lx,Ly);
v = v0*ones(Lx,Ly);

for s = 1:nsrc
    ix = round(sx0(s)/dx) + 1;
    iy = round(sy0(s)/dy) + 1;

    Cen(ix,iy) = Cen(ix,iy) + smass(s)/(dx*dy);
end

[ex,ey] = setup(e);

feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
f = feq;

rmse_v = zeros(Ns,1);
mae_v = zeros(Ns,1);
r2_v = zeros(Ns,1);
maxerr_v = zeros(Ns,1);

cmin_v = zeros(Ns,1);
neg_v = zeros(Ns,1);

mass_num = zeros(Ns,1);
mass_ana = zeros(Ns,1);
mass_err = zeros(Ns,1);

pk_num = zeros(Ns,1);
pk_ana = zeros(Ns,1);

pkx_num = zeros(Ns,1);
pky_num = zeros(Ns,1);

pkx_ana = zeros(Ns,1);
pky_ana = zeros(Ns,1);

cx_num = zeros(Ns,1);
cy_num = zeros(Ns,1);

cx_exp = zeros(Ns,1);
cy_exp = zeros(Ns,1);

hasNeg = false(Ns,1);

Cnum_all = zeros(Nx,Ny,Ns);
Cana_all = zeros(Nx,Ny,Ns);

src_final = zeros(Nx,Ny,nsrc);

fprintf('Multiple source experiment\n');
fprintf('%d sources, final time %.0f s\n\n',nsrc,max(tsnap));

si = 1;
timer = tic;

for iter = 1:Nt

    ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);
    ftemp = Peri_BC_inOut(1,Ly,ftemp,Lx);
    ftemp = Peri_BC_norSou(1,Lx,ftemp,Ly);

    [Cen,f] = solution(ftemp,1,Lx,1,Ly);

    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

    if si <= Ns && iter == iter_snap(si)

        tc = tsnap(si);
        Cnum = Cen(1:Lx-1,1:Ly-1);

        Cana = zeros(size(X));

        sx_now = zeros(nsrc,1);
        sy_now = zeros(nsrc,1);

        for s = 1:nsrc

            sx_now(s) = mod(sx0(s) + u0*tc,Lxm);
            sy_now(s) = mod(sy0(s) + v0*tc,Lym);

            Csrc = zeros(size(X));
            A = smass(s)/(4*pi*D*tc);

            for mx = -1:1
                for my = -1:1

                    xi = sx_now(s) + mx*Lxm;
                    yi = sy_now(s) + my*Lym;

                    Csrc = Csrc + A .* exp( ...
                        -((X-xi).^2 + (Y-yi).^2)/(4*D*tc));
                end
            end

            Cana = Cana + Csrc;

            if si == Ns
                src_final(:,:,s) = Csrc;
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

        merr = 100*abs(mnum-mass0)/mass0;

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

        ethx = 2*pi*sx_now/Lxm;

        eax = atan2( ...
            sum(smass.*sin(ethx)), ...
            sum(smass.*cos(ethx)) );

        if eax < 0
            eax = eax + 2*pi;
        end

        ex_c = Lxm*eax/(2*pi);

        ethy = 2*pi*sy_now/Lym;

        eay = atan2( ...
            sum(smass.*sin(ethy)), ...
            sum(smass.*cos(ethy)) );

        if eay < 0
            eay = eay + 2*pi;
        end

        ey_c = Lym*eay/(2*pi);

        rmse_v(si) = rmse;
        mae_v(si) = mae;
        r2_v(si) = r2;
        maxerr_v(si) = maxerr;

        cmin_v(si) = cmin;
        neg_v(si) = nneg;

        mass_num(si) = mnum;
        mass_ana(si) = mana;
        mass_err(si) = merr;

        pk_num(si) = pnum;
        pk_ana(si) = pana;

        pkx_num(si) = pnx;
        pky_num(si) = pny;

        pkx_ana(si) = pax;
        pky_ana(si) = pay;

        cx_num(si) = cx;
        cy_num(si) = cy;

        cx_exp(si) = ex_c;
        cy_exp(si) = ey_c;

        hasNeg(si) = nneg > 0;

        Cnum_all(:,:,si) = Cnum;
        Cana_all(:,:,si) = Cana;

        fprintf('t = %.0f s | RMSE = %.4e | R2 = %.8f | mass = %.6f\n', ...
            tc,rmse,r2,mnum);

        si = si + 1;
    end
end

runtime = toc(timer);

T = table( ...
    tsnap,rmse_v,mae_v,r2_v,maxerr_v,cmin_v,neg_v, ...
    mass_num,mass_ana,mass_err,pk_num,pk_ana,pkx_num,pky_num, ...
    pkx_ana,pky_ana,cx_num,cy_num,cx_exp,cy_exp,hasNeg, ...
    'VariableNames',{ ...
    'Time_s','RMSE','MAE','R2','MaxAbsoluteError','MinimumC', ...
    'NegativePointCount','NumericalMass','AnalyticalMass', ...
    'MassErrorPercent','NumericalPeak','AnalyticalPeak', ...
    'NumericalPeakX','NumericalPeakY','AnalyticalPeakX', ...
    'AnalyticalPeakY','NumericalCentroidX','NumericalCentroidY', ...
    'ExpectedCentroidX','ExpectedCentroidY','HasNegativeConcentration'});

fprintf('\n');
disp(T);
fprintf('Runtime = %.2f s\n',runtime);

resdir = fullfile(proj_dir,'06_Results','2D_Multiple_Sources');

if ~exist(resdir,'dir')
    mkdir(resdir);
end

writetable(T,fullfile(resdir,'2D_multiple_sources_summary.csv'));

Tsrc = table((1:nsrc)',sx0,sy0,smass, ...
    'VariableNames',{'SourceNumber','InitialX_m','InitialY_m','Mass'});

writetable(Tsrc, ...
    fullfile(resdir,'2D_multiple_sources_source_locations.csv'));

save(fullfile(resdir,'2D_multiple_sources_results.mat'));

figure;
tiledlayout(1,3);

for k = 1:Ns

    nexttile;

    contourf(xu,yu,Cnum_all(:,:,k)',20);
    hold on;

    contour(xu,yu,Cana_all(:,:,k)',8,'k--','LineWidth',1);

    xs = mod(sx0 + u0*tsnap(k),Lxm);
    ys = mod(sy0 + v0*tsnap(k),Lym);

    plot(xs,ys,'kx','MarkerSize',8,'LineWidth',1.5);

    colorbar;
    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('t = %.0f s',tsnap(k)));

    axis equal;
    xlim([0 120]);
    ylim([0 100]);

    hold off;
end

sgtitle('Transport and Interaction of Three Pollutant Sources');

exportgraphics(gcf, ...
    fullfile(resdir,'2D_multiple_sources_multitime.png'), ...
    'Resolution',300);

figure;
tiledlayout(1,2);

nexttile;

contourf(xu,yu,Cnum_all(:,:,end)',20);
colorbar;

xlabel('x (m)');
ylabel('y (m)');
title('LBM Numerical');

axis equal;
xlim([0 120]);
ylim([0 100]);

nexttile;

contourf(xu,yu,Cana_all(:,:,end)',20);
colorbar;

xlabel('x (m)');
ylabel('y (m)');
title('Analytical Superposition');

axis equal;
xlim([0 120]);
ylim([0 100]);

sgtitle(sprintf('Multiple Sources at t = %.0f s',max(tsnap)));

exportgraphics(gcf, ...
    fullfile(resdir,'2D_multiple_sources_numerical_vs_analytical.png'), ...
    'Resolution',300);

final_err = abs(Cnum_all(:,:,end)-Cana_all(:,:,end));

figure;

contourf(xu,yu,final_err',20);
colorbar;

xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Absolute Error at t = %.0f s',max(tsnap)));

axis equal;
xlim([0 120]);
ylim([0 100]);

exportgraphics(gcf, ...
    fullfile(resdir,'2D_multiple_sources_absolute_error.png'), ...
    'Resolution',300);

figure;
tiledlayout(2,2);

for s = 1:nsrc

    nexttile;

    contourf(xu,yu,src_final(:,:,s)',20);
    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('Source %d contribution',s));

    axis equal;
    xlim([0 120]);
    ylim([0 100]);
end

nexttile;

contourf(xu,yu,Cana_all(:,:,end)',20);
colorbar;

xlabel('x (m)');
ylabel('y (m)');
title('Combined concentration');

axis equal;
xlim([0 120]);
ylim([0 100]);

sgtitle('Principle of Superposition for Multiple Pollutant Sources');

exportgraphics(gcf, ...
    fullfile(resdir,'2D_multiple_sources_superposition.png'), ...
    'Resolution',300);

fprintf('Done, results in %s\n',resdir);
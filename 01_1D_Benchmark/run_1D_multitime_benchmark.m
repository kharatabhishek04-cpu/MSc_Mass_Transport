% Multi-time 1D LBM benchmark vs analytical point-source solution
% Snapshots at 5k, 10k, 15k, 20k, 25k, 30k seconds
clear; clc; close all;

% Domain
Lch = 400;
Lx = 401; Ly = 5;
dx = Lch/(Lx-1); dy = dx;
dt = 0.1;
e = dx/dt;

tau = 1.1;

u0 = 0.01; v0 = 0.0;
D = 0.01;

C0 = 1.0;
x0 = 10.0;
src_idx = round(x0/dx) + 1;

% snapshot schedule
tsnap = [5000; 10000; 15000; 20000; 25000; 30000];
iter_snap = round(tsnap/dt);
Nt = round(max(tsnap)/dt);
Nsnap = length(tsnap);

% Initialise
solid = zeros(Lx,Ly);
Cen = zeros(Lx,Ly);
u = u0*ones(Lx,Ly); v = v0*ones(Lx,Ly);
Cen(src_idx,:) = C0;

[ex,ey] = setup(e);
feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
f = feq;

xc = (0:Lx-1)'*dx;
mid_y = ceil(Ly/2);

% Pre-allocate
Cnum_all = zeros(Lx, Nsnap);
Cana_all = zeros(Lx, Nsnap);
err_all  = zeros(Lx, Nsnap);
rmse_v = zeros(Nsnap,1);  mae_v = zeros(Nsnap,1);
r2_v   = zeros(Nsnap,1);  maxerr_v = zeros(Nsnap,1);
Cmin_v = zeros(Nsnap,1);
pk_num = zeros(Nsnap,1);  pk_ana = zeros(Nsnap,1);
pkx_num = zeros(Nsnap,1); pkx_ana = zeros(Nsnap,1);
mass_num = zeros(Nsnap,1); mass_ana = zeros(Nsnap,1);

fprintf('Lx=%d, dx=%.4f, dt=%.4f, tau=%.2f\n', Lx, dx, dt, tau);
fprintf('u0=%.4f, D=%.4f, source at x=%.1f\n', u0, D, x0);
fprintf('Running %d iterations to t=%.0f s\n\n', Nt, max(tsnap));

% Main loop
si = 1;
for iter = 1:Nt

    ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);
    ftemp = Peri_BC_inOut(2,Ly,ftemp,Lx);
    ftemp = BC_NorthSouth(ftemp,1,Lx,1,Ly);

    % y-boundary copy
    for x = 1:Lx
        for a = 1:5, ftemp(a,x,1) = ftemp(a,x,2); end
    end
    for x = 1:Lx
        for a = 1:5, ftemp(a,x,Ly) = ftemp(a,x,Ly-1); end
    end

    [Cen,f] = solution(ftemp,1,Lx,1,Ly);
    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

    %fprintf('iter %d / %d\n', iter, Nt);

    % Record snapshot
    if si <= Nsnap && iter == iter_snap(si)
        tc = tsnap(si);
        Cnum = Cen(:,mid_y);
        Cana = C0./sqrt(4*pi*D*tc) .* exp(-((xc - x0 - u0*tc).^2)./(4*D*tc));
        err = Cnum - Cana;

        rmse_v(si)   = sqrt(mean(err.^2));
        mae_v(si)    = mean(abs(err));
        SS_res = sum(err.^2);
        SS_tot = sum((Cana - mean(Cana)).^2);
        r2_v(si)     = 1 - SS_res/SS_tot;
        maxerr_v(si) = max(abs(err));
        Cmin_v(si)   = min(Cnum);

        [pk_num(si), ix] = max(Cnum);  pkx_num(si) = xc(ix);
        [pk_ana(si), ix] = max(Cana);  pkx_ana(si) = xc(ix);
        mass_num(si) = sum(Cnum)*dx;
        mass_ana(si) = sum(Cana)*dx;

        Cnum_all(:,si) = Cnum;
        Cana_all(:,si) = Cana;
        err_all(:,si)  = err;

        fprintf('t=%5.0fs | RMSE=%.4e  R2=%.8f  peak: num=%.6f ana=%.6f\n', ...
            tc, rmse_v(si), r2_v(si), pk_num(si), pk_ana(si));

        si = si + 1;
    end
end

% Summary table
hasNeg = Cmin_v < -1e-12;
T = table(tsnap, rmse_v, mae_v, r2_v, maxerr_v, Cmin_v, ...
    pk_num, pk_ana, pkx_num, pkx_ana, mass_num, mass_ana, hasNeg, ...
    'VariableNames', {'Time_s','RMSE','MAE','R2','MaxAbsErr', ...
    'MinC','NumPeak','AnaPeak','NumPeakX_m','AnaPeakX_m', ...
    'NumMass','AnaMass','HasNegC'});

fprintf('\n'); disp(T);

% Save results
resdir = fullfile('..','06_Results','1D_Benchmark');
if ~exist(resdir,'dir'), mkdir(resdir); end

writetable(T, fullfile(resdir,'1D_benchmark_summary.csv'));

% profiles export (long form for plotting in Python/Excel later)
Npts = Lx * Nsnap;
prof_t = zeros(Npts,1); prof_x = zeros(Npts,1);
prof_num = zeros(Npts,1); prof_ana = zeros(Npts,1); prof_err = zeros(Npts,1);
for k = 1:Nsnap
    idx = (k-1)*Lx + (1:Lx);
    prof_t(idx) = tsnap(k);
    prof_x(idx) = xc;
    prof_num(idx) = Cnum_all(:,k);
    prof_ana(idx) = Cana_all(:,k);
    prof_err(idx) = err_all(:,k);
end
Tprof = table(prof_t, prof_x, prof_num, prof_ana, prof_err, ...
    'VariableNames', {'Time_s','x_m','C_num','C_ana','Error'});
writetable(Tprof, fullfile(resdir,'1D_benchmark_profiles.csv'));

save(fullfile(resdir,'1D_benchmark_multitime.mat'), ...
    'tsnap','xc','Cnum_all','Cana_all','err_all','T');

% Plot
figure; tiledlayout(2,3);
for k = 1:Nsnap
    nexttile;
    plot(xc, Cnum_all(:,k), 'o', 'MarkerSize', 3); hold on;
    plot(xc, Cana_all(:,k), '-', 'LineWidth', 1.5);

    xmid = x0 + u0*tsnap(k);
    xlim([max(0, xmid-80) min(Lch, xmid+80)]);
    xlabel('x (m)'); ylabel('C (kg/m)');
    title(sprintf('t = %.0f s', tsnap(k)));
    grid on;
    if k == 1, legend('LBM','Analytical'); end
    hold off;
end
sgtitle('1D LBM Multi-Time Benchmark');

%exportgraphics(gcf, fullfile(resdir,'1D_benchmark_multitime.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(resdir,'1D_benchmark_multitime.png'), 'Resolution', 300);
fprintf('Done, results in %s\n', resdir);
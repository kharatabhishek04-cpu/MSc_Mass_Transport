% 1D diffusion experiment - vary D, fixed velocity
% D = [0.006, 0.008, 0.010, 0.015, 0.020] m^2/s
% Observe at t = 10000, 20000, 25000 s
clear; clc; close all;

% Load LBM functions from benchmark folder
proj_dir = fileparts(fileparts(mfilename('fullpath')));
bm_dir = fullfile(proj_dir, '01_1D_Benchmark');
addpath(bm_dir);
fprintf('LBM functions from: %s\n\n', bm_dir);

% Domain
Lch = 400;
Lx = 401; Ly = 5;
dx = Lch/(Lx-1); dy = dx;
dt = 0.1;
e = dx/dt;

tau = 1.1;
u0 = 0.01; v0 = 0.0;

C0 = 1.0;
x0 = 10.0;
src_idx = round(x0/dx) + 1;

% experiment parameters
D_vals = [0.006; 0.008; 0.010; 0.015; 0.020];
tsnap = [10000; 20000; 25000];
iter_snap = round(tsnap/dt);
Nt = round(max(tsnap)/dt);
Nd = length(D_vals);
Ns = length(tsnap);

xc = (0:Lx-1)'*dx;
mid_y = ceil(Ly/2);

% Pre-allocate summary storage
Ncases = Nd * Ns;
res = struct();
res.caseID = strings(Ncases,1);
res.D = zeros(Ncases,1);       res.Pe = zeros(Ncases,1);
res.time = zeros(Ncases,1);
res.rmse = zeros(Ncases,1);    res.mae = zeros(Ncases,1);
res.r2 = zeros(Ncases,1);      res.maxerr = zeros(Ncases,1);
res.Cmin = zeros(Ncases,1);
res.pk_num = zeros(Ncases,1);  res.pk_ana = zeros(Ncases,1);
res.pkx_num = zeros(Ncases,1); res.pkx_ana = zeros(Ncases,1);
res.mass_num = zeros(Ncases,1); res.mass_ana = zeros(Ncases,1);
res.centroid_num = zeros(Ncases,1); res.centroid_exp = zeros(Ncases,1);
res.sigma_num = zeros(Ncases,1); res.sigma_theo = zeros(Ncases,1);
res.hasNeg = false(Ncases,1);
res.runtime = zeros(Ncases,1);

% Profile storage for plotting
Cnum_all = zeros(Lx, Nd, Ns);
Cana_all = zeros(Lx, Nd, Ns);

% long-form profile arrays
prof_D = []; prof_Pe = []; prof_t = [];
prof_x = []; prof_num = []; prof_ana = []; prof_err = [];

fprintf('Lx=%d, dx=%.4f, dt=%.4f, tau=%.2f, u0=%.4f\n', Lx, dx, dt, tau, u0);
fprintf('D values: '); fprintf('%.3f ', D_vals); fprintf('\n');
fprintf('Snapshots: '); fprintf('%.0f ', tsnap); fprintf('s\n\n');

row = 0;

for di = 1:Nd
    D = D_vals(di);
    Pe = u0*dx/D;
    fprintf('--- D=%.4f (Pe=%.3f) ---\n', D, Pe);

    % Reset model
    solid = zeros(Lx,Ly);
    Cen = zeros(Lx,Ly);
    u = u0*ones(Lx,Ly); v = v0*ones(Lx,Ly);
    Cen(src_idx,:) = C0;

    [ex,ey] = setup(e);
    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
    f = feq;

    si = 1;
    tstart = tic;

    for iter = 1:Nt
        ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);
        ftemp = Peri_BC_inOut(2,Ly,ftemp,Lx);
        ftemp = BC_NorthSouth(ftemp,1,Lx,1,Ly);

        for x = 1:Lx
            for a = 1:5, ftemp(a,x,1) = ftemp(a,x,2); end
        end
        for x = 1:Lx
            for a = 1:5, ftemp(a,x,Ly) = ftemp(a,x,Ly-1); end
        end

        [Cen,f] = solution(ftemp,1,Lx,1,Ly);
        feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

        %fprintf('D=%.3f iter %d/%d\n', D, iter, Nt);

        if si <= Ns && iter == iter_snap(si)
            tc = tsnap(si);
            row = row + 1;

            Cnum = Cen(:,mid_y);
            Cana = C0./sqrt(4*pi*D*tc) .* exp(-((xc - x0 - u0*tc).^2)./(4*D*tc));
            err = Cnum - Cana;

            % error stats
            cur_rmse = sqrt(mean(err.^2));
            cur_mae  = mean(abs(err));
            SS_res = sum(err.^2);
            SS_tot = sum((Cana - mean(Cana)).^2);
            cur_r2 = 1 - SS_res/SS_tot;
            cur_maxerr = max(abs(err));
            cur_Cmin = min(Cnum);

            [cur_pk_num, ix] = max(Cnum); cur_pkx_num = xc(ix);
            [cur_pk_ana, ix] = max(Cana); cur_pkx_ana = xc(ix);
            cur_mass_num = sum(Cnum)*dx;
            cur_mass_ana = sum(Cana)*dx;

            % centroid and spread
            cur_centroid = sum(xc .* Cnum)*dx / cur_mass_num;
            exp_centroid = x0 + u0*tc;
            cur_sigma = sqrt(sum(((xc - cur_centroid).^2) .* Cnum)*dx / cur_mass_num);
            theo_sigma = sqrt(2*D*tc);

            neg_flag = cur_Cmin < -1e-12;
            cur_runtime = toc(tstart);

            % Store in struct
            res.caseID(row) = sprintf('D%.3f_T%d', D, tc);
            res.D(row) = D; res.Pe(row) = Pe; res.time(row) = tc;
            res.rmse(row) = cur_rmse; res.mae(row) = cur_mae;
            res.r2(row) = cur_r2; res.maxerr(row) = cur_maxerr;
            res.Cmin(row) = cur_Cmin;
            res.pk_num(row) = cur_pk_num; res.pk_ana(row) = cur_pk_ana;
            res.pkx_num(row) = cur_pkx_num; res.pkx_ana(row) = cur_pkx_ana;
            res.mass_num(row) = cur_mass_num; res.mass_ana(row) = cur_mass_ana;
            res.centroid_num(row) = cur_centroid; res.centroid_exp(row) = exp_centroid;
            res.sigma_num(row) = cur_sigma; res.sigma_theo(row) = theo_sigma;
            res.hasNeg(row) = neg_flag;
            res.runtime(row) = cur_runtime;

            Cnum_all(:,di,si) = Cnum;
            Cana_all(:,di,si) = Cana;

            % append to long-form profiles
            prof_D   = [prof_D;   repmat(D,Lx,1)];
            prof_Pe  = [prof_Pe;  repmat(Pe,Lx,1)];
            prof_t   = [prof_t;   repmat(tc,Lx,1)];
            prof_x   = [prof_x;   xc];
            prof_num = [prof_num; Cnum];
            prof_ana = [prof_ana; Cana];
            prof_err = [prof_err; err];

            fprintf('  t=%5.0fs | RMSE=%.4e  R2=%.8f  sigma: num=%.3f theo=%.3f\n', ...
                tc, cur_rmse, cur_r2, cur_sigma, theo_sigma);

            si = si + 1;
        end
    end
end

% Build summary table
T = table(res.caseID, res.D, res.Pe, res.time, res.rmse, res.mae, res.r2, ...
    res.maxerr, res.Cmin, res.pk_num, res.pk_ana, res.pkx_num, res.pkx_ana, ...
    res.mass_num, res.mass_ana, res.centroid_num, res.centroid_exp, ...
    res.sigma_num, res.sigma_theo, res.hasNeg, res.runtime, ...
    'VariableNames', {'CaseID','Diffusion','Peclet','Time_s','RMSE','MAE','R2', ...
    'MaxAbsErr','MinC','NumPeak','AnaPeak','NumPeakX_m','AnaPeakX_m', ...
    'NumMass','AnaMass','NumCentroidX_m','ExpCentroidX_m', ...
    'NumSigma_m','TheoSigma_m','HasNegC','Runtime_s'});

fprintf('\n'); disp(T);

% save results
resdir = fullfile(proj_dir, '06_Results', '1D_Diffusion');
if ~exist(resdir,'dir'), mkdir(resdir); end

writetable(T, fullfile(resdir, '1D_diffusion_summary.csv'));

Tprof = table(prof_D, prof_Pe, prof_t, prof_x, prof_num, prof_ana, prof_err, ...
    'VariableNames', {'Diffusion','Peclet','Time_s','x_m','C_num','C_ana','Error'});
writetable(Tprof, fullfile(resdir, '1D_diffusion_profiles.csv'));

save(fullfile(resdir, '1D_diffusion_experiment.mat'));

% Plot: effect of D at each snapshot time
figure; tiledlayout(1,3);
for ti = 1:Ns
    nexttile; hold on;
    for di = 1:Nd
        plot(xc, Cnum_all(:,di,ti), 'LineWidth', 1.5);
    end
    xmid = x0 + u0*tsnap(ti);
    xlim([max(0, xmid-120) min(Lch, xmid+120)]);
    xlabel('x (m)'); ylabel('C (kg/m)');
    title(sprintf('t = %.0f s', tsnap(ti)));
    grid on;
    if ti == 1, legend('D=0.006','D=0.008','D=0.010','D=0.015','D=0.020'); end
    hold off;
end
sgtitle('1D Controlled Diffusion Experiment');
exportgraphics(gcf, fullfile(resdir, '1D_diffusion_comparison.png'), 'Resolution', 300);

% Plot: num vs analytical at t=20000s
val_ti = 2;  % index for t=20000
figure; tiledlayout(2,3);
for di = 1:Nd
    nexttile;
    plot(xc, Cnum_all(:,di,val_ti), 'o', 'MarkerSize', 3); hold on;
    plot(xc, Cana_all(:,di,val_ti), '-', 'LineWidth', 1.5);
    xmid = x0 + u0*tsnap(val_ti);
    xlim([max(0, xmid-100) min(Lch, xmid+100)]);
    xlabel('x (m)'); ylabel('C (kg/m)');
    title(sprintf('D = %.3f', D_vals(di)));
    grid on;
    if di == 1, legend('LBM','Analytical'); end
    hold off;
end
sgtitle('Numerical vs Analytical at t = 20000 s');
%exportgraphics(gcf, fullfile(resdir, '1D_diffusion_validation_20000s.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(resdir, '1D_diffusion_validation_20000s.png'), 'Resolution', 300);

fprintf('Done, results in %s\n', resdir);
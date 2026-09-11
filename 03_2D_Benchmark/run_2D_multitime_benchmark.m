clear; clc; close all;

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
D = 0.02;

C0 = 1.0;
x0 = 10.0;
y0 = 50.0;

ix0 = round(x0/dx) + 1;
iy0 = round(y0/dy) + 1;

tsnap = [2000; 4000; 6000; 8000; 10000; 12000];
iter_snap = round(tsnap/dt);
Ns = length(tsnap);
Nt = round(max(tsnap)/dt);

solid = zeros(Lx,Ly);
Cen = zeros(Lx,Ly);
u = u0*ones(Lx,Ly);
v = v0*ones(Lx,Ly);

Cen(ix0,iy0) = C0;

[ex,ey] = setup(e);
feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
f = feq;

xu = (0:Lx-2)'*dx;
yu = (0:Ly-2)*dy;
[X,Y] = ndgrid(xu,yu);

Nx = length(xu);
Ny = length(yu);

rmse_v = zeros(Ns,1);
mae_v = zeros(Ns,1);
r2_v = zeros(Ns,1);
maxerr_v = zeros(Ns,1);
minC_v = zeros(Ns,1);
neg_v = zeros(Ns,1);

mass_num = zeros(Ns,1);
mass_ana = zeros(Ns,1);

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

sigx = zeros(Ns,1);
sigy = zeros(Ns,1);
sig_theory = zeros(Ns,1);

hasNeg = false(Ns,1);

Cnum_all = zeros(Nx,Ny,Ns);
Cana_all = zeros(Nx,Ny,Ns);

npts = Nx*Ny;
nrows = npts*Ns;

prof_t = zeros(nrows,1);
prof_x = zeros(nrows,1);
prof_y = zeros(nrows,1);
prof_num = zeros(nrows,1);
prof_ana = zeros(nrows,1);
prof_err = zeros(nrows,1);

p0 = 1;

fprintf('2D benchmark\n');
fprintf('Domain = %.0f x %.0f m\n',Lxm,Lym);
fprintf('Grid = %d x %d\n',Lx,Ly);
fprintf('Validation grid = %d x %d\n',Nx,Ny);
fprintf('dx = %.4f, dy = %.4f, dt = %.4f\n',dx,dy,dt);
fprintf('u = %.4f, v = %.4f, D = %.4f, tau = %.2f\n',u0,v0,D,tau);
fprintf('Source = (%.2f, %.2f) m\n',x0,y0);
fprintf('Running %d iterations\n\n',Nt);

si = 1;

for iter = 1:Nt

    ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);
    ftemp = Peri_BC_inOut(1,Ly,ftemp,Lx);
    ftemp = Peri_BC_norSou(1,Lx,ftemp,Ly);

    [Cen,f] = solution(ftemp,1,Lx,1,Ly);
    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

    if si <= Ns && iter == iter_snap(si)

        tc = tsnap(si);
        Cnum = Cen(1:Lx-1,1:Ly-1);

        xexp = mod(x0 + u0*tc,Lxm);
        yexp = mod(y0 + v0*tc,Lym);

        Cana = zeros(size(X));
        A = C0/(4*pi*D*tc);

        for mx = -1:1
            for my = -1:1
                xi = xexp + mx*Lxm;
                yi = yexp + my*Lym;

                Cana = Cana + A .* exp( ...
                    -((X-xi).^2 + (Y-yi).^2)/(4*D*tc));
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

        [pnum,ind] = max(Cnum(:));
        [inx,iny] = ind2sub(size(Cnum),ind);
        pnx = xu(inx);
        pny = yu(iny);

        [pana,ind] = max(Cana(:));
        [iax,iay] = ind2sub(size(Cana),ind);
        pax = xu(iax);
        pay = yu(iay);

        thx = 2*pi*X/Lxm;
        ax = atan2(sum(Cnum(:).*sin(thx(:))), ...
            sum(Cnum(:).*cos(thx(:))));

        if ax < 0
            ax = ax + 2*pi;
        end

        cx = Lxm*ax/(2*pi);

        thy = 2*pi*Y/Lym;
        ay = atan2(sum(Cnum(:).*sin(thy(:))), ...
            sum(Cnum(:).*cos(thy(:))));

        if ay < 0
            ay = ay + 2*pi;
        end

        cy = Lym*ay/(2*pi);

        distx = mod(X-cx+Lxm/2,Lxm)-Lxm/2;
        disty = mod(Y-cy+Lym/2,Lym)-Lym/2;

        sx = sqrt(sum((distx(:).^2).*Cnum(:))*dx*dy/mnum);
        sy = sqrt(sum((disty(:).^2).*Cnum(:))*dx*dy/mnum);
        st = sqrt(2*D*tc);

        rmse_v(si) = rmse;
        mae_v(si) = mae;
        r2_v(si) = r2;
        maxerr_v(si) = maxerr;
        minC_v(si) = cmin;
        neg_v(si) = nneg;

        mass_num(si) = mnum;
        mass_ana(si) = mana;

        pk_num(si) = pnum;
        pk_ana(si) = pana;

        pkx_num(si) = pnx;
        pky_num(si) = pny;
        pkx_ana(si) = pax;
        pky_ana(si) = pay;

        cx_num(si) = cx;
        cy_num(si) = cy;
        cx_exp(si) = xexp;
        cy_exp(si) = yexp;

        sigx(si) = sx;
        sigy(si) = sy;
        sig_theory(si) = st;

        hasNeg(si) = nneg > 0;

        Cnum_all(:,:,si) = Cnum;
        Cana_all(:,:,si) = Cana;

        p1 = p0 + npts - 1;

        prof_t(p0:p1) = tc;
        prof_x(p0:p1) = X(:);
        prof_y(p0:p1) = Y(:);
        prof_num(p0:p1) = Cnum(:);
        prof_ana(p0:p1) = Cana(:);
        prof_err(p0:p1) = err(:);

        p0 = p1 + 1;

        fprintf('t = %.0f s\n',tc);
        fprintf('RMSE = %.12e, MAE = %.12e, R2 = %.10f\n',rmse,mae,r2);
        fprintf('Max error = %.12e, Min C = %.12e, negative = %d\n', ...
            maxerr,cmin,nneg);
        fprintf('Mass: num = %.10f, ana = %.10f\n',mnum,mana);
        fprintf('Peak: num = %.10f at (%.2f, %.2f), ana = %.10f at (%.2f, %.2f)\n', ...
            pnum,pnx,pny,pana,pax,pay);
        fprintf('Centre: num = (%.4f, %.4f), expected = (%.4f, %.4f)\n', ...
            cx,cy,xexp,yexp);
        fprintf('Sigma: x = %.6f, y = %.6f, theory = %.6f\n\n',sx,sy,st);

        si = si + 1;
    end
end

T = table(tsnap,rmse_v,mae_v,r2_v,maxerr_v,minC_v,neg_v, ...
    mass_num,mass_ana,pk_num,pk_ana,pkx_num,pky_num,pkx_ana,pky_ana, ...
    cx_num,cy_num,cx_exp,cy_exp,sigx,sigy,sig_theory,hasNeg, ...
    'VariableNames',{'Time_s','RMSE','MAE','R2','MaxAbsoluteError', ...
    'MinimumC','NegativePointCount','NumericalMass','AnalyticalMass', ...
    'NumericalPeak','AnalyticalPeak','NumericalPeakX','NumericalPeakY', ...
    'AnalyticalPeakX','AnalyticalPeakY','NumericalCentroidX', ...
    'NumericalCentroidY','ExpectedCentroidX','ExpectedCentroidY', ...
    'SigmaX','SigmaY','TheoreticalSigma','HasNegativeConcentration'});

fprintf('\n');
disp(T);

project_folder = fileparts(fileparts(mfilename('fullpath')));
resdir = fullfile(project_folder,'06_Results','2D_Benchmark');

if ~exist(resdir,'dir')
    mkdir(resdir);
end

writetable(T,fullfile(resdir,'2D_benchmark_summary.csv'));

Tprof = table(prof_t,prof_x,prof_y,prof_num,prof_ana,prof_err, ...
    'VariableNames',{'Time_s','x_m','y_m','NumericalConcentration', ...
    'AnalyticalConcentration','Error'});

writetable(Tprof,fullfile(resdir,'2D_benchmark_profiles.csv'));

save(fullfile(resdir,'2D_benchmark_multitime.mat'));

figure;
tiledlayout(2,3);

for k = 1:Ns
    nexttile;

    contourf(xu,yu,Cnum_all(:,:,k)',20);
    hold on;
    contour(xu,yu,Cana_all(:,:,k)',8,'k--','LineWidth',1);

    colorbar;
    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('t = %.0f s',tsnap(k)));

    axis equal;
    xlim([0 Lxm]);
    ylim([0 Lym]);
    hold off;
end

sgtitle('2D LBM Numerical Field with Analytical Contours');

exportgraphics(gcf,fullfile(resdir,'2D_benchmark_multitime_contours.png'), ...
    'Resolution',300);

[~,iy_mid] = min(abs(yu-y0));

figure;
tiledlayout(2,3);

for k = 1:Ns
    nexttile;

    plot(xu,Cnum_all(:,iy_mid,k),'o','MarkerSize',3);
    hold on;
    plot(xu,Cana_all(:,iy_mid,k),'-','LineWidth',1.5);

    xlabel('x (m)');
    ylabel('C');
    title(sprintf('t = %.0f s',tsnap(k)));
    grid on;

    if k == 1
        legend('LBM Numerical','Analytical');
    end

    hold off;
end

sgtitle('2D Benchmark Centreline Validation at y = 50 m');

exportgraphics(gcf,fullfile(resdir,'2D_benchmark_centreline_validation.png'), ...
    'Resolution',300);

fprintf('\nDone, results in %s\n',resdir);
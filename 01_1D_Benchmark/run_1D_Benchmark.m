clear; clc; close all;

% Domain
Lch = 400;                 % channel length (m)
Lx = 401; Ly = 5;
dx = Lch/(Lx-1); dy = dx;
dt = 0.1;
e = dx/dt;

tau = 1.1;

% Flow and diffusion
u0 = 0.01; v0 = 0.0;
D = 0.01;

% Source
C0 = 1.0;
x0 = 10.0;
src_idx = round(x0/dx) + 1;

tfinal = 10;               % seconds
Nt = round(tfinal/dt);

% Initialise arrays
solid = zeros(Lx,Ly);
Cen = zeros(Lx,Ly);
u = u0*ones(Lx,Ly); v = v0*ones(Lx,Ly);

Cen(src_idx,:) = C0;

[ex,ey] = setup(e);

feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
f = feq;

fprintf('Lx=%d, dx=%.4f, dt=%.4f, tau=%.2f\n', Lx, dx, dt, tau);
fprintf('u0=%.4f, D=%.4f, source at x=%.1f (idx %d)\n', u0, D, x0, src_idx);
fprintf('Running %d iterations (%.1f s)\n\n', Nt, tfinal);

% Main loop
for iter = 1:Nt

    ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);

    ftemp = Peri_BC_inOut(2,Ly,ftemp,Lx);
    ftemp = BC_NorthSouth(ftemp,1,Lx,1,Ly);

    % Copy lower/upper y-boundaries
    for x = 1:Lx
        for a = 1:5, ftemp(a,x,1) = ftemp(a,x,2); end
    end
    for x = 1:Lx
        for a = 1:5, ftemp(a,x,Ly) = ftemp(a,x,Ly-1); end
    end

    [Cen,f] = solution(ftemp,1,Lx,1,Ly);
    %fprintf('iter %d\n', iter);

    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
end

% Extract centreline profile
tsim = Nt*dt;
mid_y = ceil(Ly/2);
Cnum = Cen(:,mid_y);

fprintf('Done. t=%.2f s, Cmin=%.6e, Cmax=%.6e, mass=%.6f\n', ...
    tsim, min(Cnum), max(Cnum), sum(Cnum)*dx);

% Analytical solution (1D instantaneous point source)
xc = (0:Lx-1)' * dx;
Cana = C0 ./ sqrt(4*pi*D*tsim) .* ...
    exp(-((xc - x0 - u0*tsim).^2) ./ (4*D*tsim));

% error metrics
err = Cnum - Cana;
RMSE = sqrt(mean(err.^2));
MAE  = mean(abs(err));
SS_res = sum(err.^2);
SS_tot = sum((Cana - mean(Cana)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf('RMSE=%.6e, MAE=%.6e, max|err|=%.6e\n', RMSE, MAE, max(abs(err)));
fprintf('R2 = %f\n', R2);

% Peak comparison
[pk_num, idx_num] = max(Cnum);
[pk_ana, idx_ana] = max(Cana);
fprintf('Peak: num=%.6f @ x=%.2f m, ana=%.6f @ x=%.2f m\n', ...
    pk_num, xc(idx_num), pk_ana, xc(idx_ana));

% Plot
figure;
plot(xc, Cnum, 'o', 'MarkerSize', 4); hold on;
plot(xc, Cana, '-', 'LineWidth', 1.5);
xlabel('Distance x (m)'); ylabel('Concentration C (kg/m)');
title(sprintf('1D LBM Validation at t = %.0f s', tsim));
legend('LBM Numerical','Analytical Solution');
grid on; xlim([0 30]); hold off;
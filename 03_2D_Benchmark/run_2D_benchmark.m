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
v0 = 0.0;
D = 0.02;

C0 = 1.0;
x0 = 10.0;
y0 = 50.0;

ix0 = round(x0/dx) + 1;
iy0 = round(y0/dy) + 1;

tfinal = 10;
Nt = round(tfinal/dt);

solid = zeros(Lx,Ly);
Cen = zeros(Lx,Ly);

u = u0*ones(Lx,Ly);
v = v0*ones(Lx,Ly);

Cen(ix0,iy0) = C0;

[ex,ey] = setup(e);

feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
f = feq;

fprintf('Running 2D test for %.1f s (%d iterations)\n', tfinal,Nt);

for iter = 1:Nt

    ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);

    ftemp = Peri_BC_inOut(1,Ly,ftemp,Lx);
    ftemp = Peri_BC_norSou(1,Lx,ftemp,Ly);

    [Cen,f] = solution(ftemp,1,Lx,1,Ly);

    feq = compute_feq(u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);
end

tsim = Nt*dt;

Cmin = min(Cen(:));
Cmax = max(Cen(:));
nneg = sum(Cen(:) < -1e-12);

mass = sum(Cen(:))*dx*dy;

[pk,idx] = max(Cen(:));
[ix,iy] = ind2sub(size(Cen),idx);

pkx = (ix-1)*dx;
pky = (iy-1)*dy;

xexp = x0 + u0*tsim;
yexp = y0 + v0*tsim;

xc = (0:Lx-1)'*dx;
yc = (0:Ly-1)*dy;

[X,Y] = ndgrid(xc,yc);

cx = sum(X(:).*Cen(:))*dx*dy/mass;
cy = sum(Y(:).*Cen(:))*dx*dy/mass;

fprintf('\nt = %.2f s\n', tsim);
fprintf('Cmin = %.12e, Cmax = %.10f\n', Cmin,Cmax);
fprintf('Negative points = %d\n', nneg);
fprintf('Mass = %.10f\n', mass);
fprintf('Peak = %.10f at (%.2f, %.2f) m\n', pk,pkx,pky);
fprintf('Centroid = (%.6f, %.6f) m\n', cx,cy);
fprintf('Expected centre = (%.6f, %.6f) m\n', xexp,yexp);

figure;

contourf(xc,yc,Cen',20);
colorbar;

xlabel('x (m)');
ylabel('y (m)');
title(sprintf('2D Pollutant Concentration at t = %.0f s',tsim));

axis equal;
xlim([0 30]);
ylim([35 65]);
grid on;
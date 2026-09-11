
clear all
load result.mat

for x=1:Lx
    xn(x)=(x-1)*dx;
    yn(x)=Cen(x,1);
end
M=numel(xn);
n=1;
index=1:n:M;
plot(xn(index),yn(index),'ro')

hold on

C0=1;
u=0.01;
D=0.01;

x0=(10-1)*dx;

xa = 0:400/50:(Lx-1)*dy;

t1=10000;
t2=30000;

ya = C0/sqrt(4*3.14*D*t1)*exp( -( (xa-x0-u*t1).^2/(4*D*t1) ) );

plot(xa,ya,'b')
xlabel('x(m)'), ylabel('C(kg/m)')
title ('Concentration')
legend 'Numerical' 'Analytical'
hold off
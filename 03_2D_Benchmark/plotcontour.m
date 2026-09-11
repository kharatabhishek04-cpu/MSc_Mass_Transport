% Plot contour
clear
clc
% The numerical data
load result.mat


x=1:Lx;
y=1:Ly;

% [xx,yy]=meshgrid(y,x);
% contour(xx,yy,Cen,10);

% use a transpose to get x y corret
[xx,yy]=meshgrid(x,y);
contour(xx,yy,Cen',10);  

xlabel('x(m)'), ylabel('y(m)')
title (' Contours of the Water head.')

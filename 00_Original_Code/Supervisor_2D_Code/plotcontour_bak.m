% Plot contour

% The numerical data
load result.mat



x=1:Lx;
y=1:Ly;
% [xx,yy]=meshgrid(x,y);
[xx,yy]=meshgrid(y,x);

z = Cen(x,y);
contour(xx,yy,z,10);


xlabel('x(m)'), ylabel('y(m)')
title (' Contours of the Water head.')

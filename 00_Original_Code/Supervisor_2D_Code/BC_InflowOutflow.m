%BC_InflowOutflow (Xs,Xe,Ys,Ye)
function [ftemp]=BC_InflowOutflow(ftemp,Xs,Xe,Ys,Ye)
x = Xs;
for y = Ys:Ye
    ftemp(1,x,y) = ftemp(3,x,y);
end

x = Xe;
for y = Ys:Ye
    ftemp(3,x,y) = ftemp(1,x,y);
end
end
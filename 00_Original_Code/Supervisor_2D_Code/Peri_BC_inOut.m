%BC_periodX

function [ftemp] = Peri_BC_inOut(Ys,Ye,ftemp,Lx) %does ftemp

%need defining

for y = Ys:Ye
    ftemp(1,1,y)=ftemp(1,Lx,y);
    ftemp(3,Lx,y)=ftemp(3,1,y);
end
    


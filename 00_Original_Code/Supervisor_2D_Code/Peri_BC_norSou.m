%BC_periodX

function [ftemp] = Peri_BC_norSou(Xs,Xe,ftemp,Ly) %does ftemp

%need defining

for x = Xs:Xe
    ftemp(2,x,1)=ftemp(2,x,Ly);
    ftemp(4,x,Ly)=ftemp(4,x,1);
end
    


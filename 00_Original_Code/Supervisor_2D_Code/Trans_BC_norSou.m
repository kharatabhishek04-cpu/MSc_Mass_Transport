%BC_periodX

function [ftemp] = Trans_BC_norSou(Xs,Xe,ftemp,Ly) %does ftemp

%need defining

for x = Xs:Xe
    ftemp(2,x,1)=ftemp(2,x,2);
    ftemp(4,x,Ly)=ftemp(4,x,Ly-1);
end
    


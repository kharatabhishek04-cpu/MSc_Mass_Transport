%BC_NorthSouth (Xs,Xe,Ys,Ye)
function[ftemp]=BC_NorthSouth(ftemp,Xs,Xe,Ys,Ye)    
y = Ye; 
for x = Xs:Xe
    ftemp(4,x,y) = ftemp(2,x,y);
end
y = Ys;
for x = Xs:Xe
    ftemp(2,x,y) = ftemp(4,x,y);
end
end
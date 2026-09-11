%solution (Xs,Xe,Ys,Ye)
function[Cen,f]=solution(ftemp,Xs,Xe,Ys,Ye)
for y = Ys:Ye
    for x = Xs:Xe
        for a=1:5;
            f(a,x,y)=ftemp(a,x,y);
        end
    end
end
for y = Ys:Ye
    for x = Xs: Xe        
        Cen(x,y)=0;
        for a = 1:5
            Cen(x,y) = Cen(x,y)+f(a,x,y);
        end
    end
end
end
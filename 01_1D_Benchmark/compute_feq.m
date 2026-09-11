%Set the initial distribution function to feq
function [feq]=compute_feq(u,v,dt,tau,e,ex,ey,Dr,solid,Lx,Ly,Cen)

lamda = Dr/(2.*dt*(tau-0.5)*e*e);

ex0 = e;
ey0 = e;
 
for y = 1:Ly
    for x = 1:Lx
        if (solid(x,y) < 1)
            for a = 1:4;
                if (a == 1 || a == 3)
             feq(a,x,y) = ey0/ex0*lamda*Cen(x,y)+Cen(x,y)/(2.*ex0*ex0)*...
                        (ex(a)*u(x,y)+ey(a)*v(x,y));
                else
             feq(a,x,y) = ex0/ey0*lamda*Cen(x,y)+Cen(x,y)/(2.*ex0*ex0)*...
                        (ex(a)*u(x,y)+ey(a)*v(x,y));
                end
            end
            feq(5,x,y) = (1.-2.*(ey0/ex0+ex0/ex0)*lamda)*Cen(x,y);
        end
    end
end
end
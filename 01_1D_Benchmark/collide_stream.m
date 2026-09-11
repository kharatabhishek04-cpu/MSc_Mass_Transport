%collide_stream (Xs,Xe,Ys,Ye)
    function [ftemp]=collide_stream(Xs,Xe,Ys,Ye,f,feq,tau)
    for y = Ys:Ye
       yp = y + 1;
       yn = y - 1 ;
       for x = Xs:Xe
          xp = x + 1;
          xn = x - 1;
          if xp <= Xe;
               ftemp(1,xp,y) = f(1,x,y) - (f(1,x,y)-feq(1,x,y))/tau;
          end
          if yp <= Ye;
               ftemp(2,x,yp) = f(2,x,y) - (f(2,x,y)-feq(2,x,y))/tau;
          end
          if xn >= Xs;
               ftemp(3,xn,y) = f(3,x,y) - (f(3,x,y)-feq(3,x,y))/tau;
          end
          if yn >= Ys;
               ftemp(4,x,yn) = f(4,x,y) - (f(4,x,y)-feq(4,x,y))/tau;
          end
               ftemp(5,x,y) = f(5,x,y) - (f(5,x,y)-feq(5,x,y))/tau;
       end
    end
    end

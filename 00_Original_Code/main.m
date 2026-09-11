%main program

% global solid Dr

max_t=1;            

length = 400;

uo = 0.01;

Lx = 401; Ly = 5;

solid = zeros(Lx,Ly);
feq = zeros(5,Lx,Ly);

dx = real(length/(Lx-1));
dy = dx;  %dx = 2 dx must be bigger than sqrt(2dt/(tau-0.5))

width = (Ly-1)*dy;

dt = dx/10;

e=dx/dt;

width = (Ly-1)*dy;

tau = 1.1;


% given velocity field
for x = 1:Lx
    for y = 1:Ly
        u(x,y) = 0.01;
        v(x,y) = 0.;
    end
end


  % initialize the depth and velocity field

  Dr = 0.01;

  for y = 1:Ly
      for x = 1:Lx
          Cen(x,y)=0.;
          if (x == 10 ) Cen(x,y) = 1.;
          end
      end
  end

%call setup
[ex,ey]=setup(e);


fprintf('%f  %f  %f  %f\n', Lx, Ly, e, dt );

%compute the equilibrium distribution function feq
[feq]=compute_feq(u,v,dt,tau,e,ex,ey,Dr,solid,Lx,Ly,Cen);

%set the initial distribution function to feq
for y=1:Ly
    for x=1:Lx
        for a=1:5
            f(a,x,y)=feq(a,x,y);
        end
    end
end

time=0;

iteration=0;

%main loop for time marching

while iteration < max_t
    
    iteration=iteration+1;
    
    time=iteration*dt;
    
    
    %call collide_stream (1,Lx,1,Ly)
    [ftemp]=collide_stream(1,Lx,1,Ly,f,feq,tau);
    
    [ftemp] = Peri_BC_inOut(2,Ly,ftemp,Lx);
    
    % call BC_NorthSouth (1,Lx,1,Ly)
    [ftemp]=BC_NorthSouth(ftemp,1,Lx,1,Ly);
    
    
    
    y = 1;
    for x = 1:Lx
        for a = 1:5
            ftemp(a,x,y) = ftemp(a,x,y+1);
        end
    end
    
    y = Ly;
    for x = 1:Lx
        for a=1:5
            ftemp(a,x,y) = ftemp(a,x,y-1);
        end
    end
    
    [Cen,f]=solution(ftemp,1,Lx,1,Ly);
    
%     x = 1;
%     for y = 1:Ly
%         for a = 1:5
%             f(a,x,y) = f(a,x+1,y);
%         end
%     end
%     
%     x = Lx;
%     for y = 1:Ly
%         for a = 1:5
%             f(a,x,y) = f(a,x-1,y);
%         end
%     end
    
   
    %update the feq
     [feq]=compute_feq(u,v,dt,tau,e,ex,ey,Dr,solid,Lx,Ly,Cen);
     
    
    % disp(sprintf('%d  %f %f', iteration, Cen(int16(Lx/2),int16(Ly/2)), time ));
    fprintf('%d  %f %f\n', iteration, Cen(int16(Lx/2),int16(Ly/2)), time );
    
    if iteration == max_t
        
        %write file
        Mat=write_file(Lx,Ly,dx,dy,dt,Cen,u,v,tau,iteration,length,width,max_t);

        save('result');
        
        count = input('Contiue computations (y/n)? ', 's');
        
        if count == 'y';
            addItea = input('Type additional number: ');
            max_t = max_t+addItea;
        end
    end

end


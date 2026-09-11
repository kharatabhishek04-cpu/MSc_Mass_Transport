function Mat=write_file(Lx,Ly,dx,dy,dt,Cen,u,v,tau,iteration,length,width,max_t)

W=fopen('result.dat','w+');

fprintf(W,' # %s',datestr(now));
fprintf(W,'\n # ......................................................................................');
fprintf(W,'\n #  length=');
fprintf(W,' %10.5f',length);
fprintf(W,' width=');
fprintf(W,' %10.5f',width);
fprintf(W,'\n #  Lx=');
fprintf(W,' %d',Lx);
fprintf(W,' Ly=');
fprintf(W,' %d',Ly);
fprintf(W,' dx=');
fprintf(W,' %10.5f',dx);
fprintf(W,' dy=');
fprintf(W,' %10.5f',dy);

fprintf(W,'\n #  dt=');
fprintf(W,' %10.5f',dt);
fprintf(W,' Total iterations = ');
fprintf(W,'%d\n',iteration);

fprintf(W,'\n # tau=');
fprintf(W,'%10.5f',tau);

fprintf(W,'\n # ...............................\n');
fprintf(W,' # Result of Computations');
fprintf(W,'\n # ...............................\n');

Ch=['    x', '          y', '           Cen', '           u', '           v'];
fprintf(W,' # %13s  %16s  %19s  %16s  %19s\n',Ch);

for j=1:Ly
    for i=1:Lx;
%        for j=1:Ly
        a=Cen(i,j);
        x=(i-1)*dx;
        y=(j-1)*dy; 
        Mat=[x,y,a,u(i,j),v(i,j)];
        fprintf(W,'\n%10.5f  %10.5f  %10.6e %10.5f  %10.5f',Mat);
    end    
    fprintf(W,'\n'); 
end
fclose(W);
return
        
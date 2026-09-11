set terminal windows
set pointsize 1.
set size 0.96,0.6
set size ratio 0.4
set key 350,0.036

ho=1.
stepy=1
#stepy=0.005
#stepy=0.0025
#stepy=0.002
#stepy=0.00125

xmin=0
xmax=400

ymin=0.
ymax=0.04
#ymax=0.1

dx=1
x0 = (10-2)*dx
Dr = 0.01
  

set xtic xmin,50,xmax
set ytic ymin,0.01,ymax
set mxtics 5
set mytics 5
#set grid mxtics mytics
#set title "This is circle segment"
set xlabel "x (m)"
set ylabel "C (kg/m)"

set nolabel

set label 't=10000s' at 80,0.032
set label 't=30000s' at 280,0.02


C0=1
u=0.01
D=0.01

x0=(10-1)*dx

t1=10000
t2=30000

sol1(x)=C0/sqrt(4*3.14*D*t1)*exp( -( (x-x0-u*t1)**2/(4*D*t1) ) )

sol2(x)=C0/sqrt(4*3.14*D*t2)*exp( -( (x-x0-u*t2)**2/(4*D*t2) ) )



plot [xmin:xmax] [ymin:ymax] "result-t10000.dat" u (int($2/stepy)==2? ($1-0):1/0):3 t "Numerical" w l 3,\
     sol1(x) t "Exact" w l 1,\
     "result-t30000.dat" u (int($2/stepy)==2? ($1-0):1/0):3 t "" w l 3,\
     sol2(x) t ""  w l 1


	
set terminal postscript eps enhanced 20
set output "com.ps"
rep
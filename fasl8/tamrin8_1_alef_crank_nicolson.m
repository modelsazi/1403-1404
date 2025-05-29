clc;clear;close all

Ti = 300; % initial temperature
p = 960; Cp = 2200; R = 0.3; h = 75; qdot = 1000; Tinf = 300;
k = 5;

ndx = 40; nx = ndx+1;
dr = R/ndx;
r = 0:dr:R;

ndt=1200; Lt = 4000;
nt=ndt+1;
dt=Lt/ndt;
t=0:dt:Lt;

alfa = k*dt/2/p/Cp/dr^2;
beta = qdot*dt/p/Cp;

T = ones(nt,nx)*Ti;

for m = 1:nt-1
    A = zeros(nx,nx); b = zeros(nx,1);
    A(1,1) = 1; A(1,2) = -1; b(1) = 0;
    for i = 2:nx-1
        A(i,i-1) = -alfa;
        A(i,i)   = 1 + alfa/r(i)*( r(i+1) + r(i) );
        A(i,i+1) = -alfa/r(i)*r(i+1);
        b(i) = alfa*T(m,i-1)+( 1 - alfa/r(i)*( r(i+1) + r(i)) )*T(m,i)+alfa/r(i)*r(i+1)*T(m,i+1)+beta;
    end
    A(nx,nx) = 1 + h*dr/k; A(nx,nx-1) = -1; b(nx) = h*dr/k*Tinf;
    T(m+1,:) = (A\b)';
end
mesh(r,t,T)
xlabel('r (m)')
ylabel('t (s)')
zlabel('T (^oC)')


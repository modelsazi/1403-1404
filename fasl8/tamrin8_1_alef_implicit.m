clc;clear;close all

Ti = 300; % initial temperature
p = 960; Cp = 2200; R = 0.3; h = 75; qdot = 1000; Tinf = 300;
k = 5;

ndx = 40; nx = ndx+1;
dr = R/ndx;
r = 0:dr:R;

ndt=100; Lt = 4000;
nt=ndt+1;
dt=Lt/ndt;
t=0:dt:Lt;

alpha=k/p/Cp;

T=zeros(nt,nx);

T(1,:)=Ti;

for m=1:nt-1
    A=zeros(nx,nx);
    b=zeros(nx,1);
    for i=2:nx-1
        A(i,i)=-1/dt-2*alpha/dr^2-alpha/r(i)/dr;
        A(i,i+1)= alpha/dr^2 + alpha/dr/r(i);
        A(i,i-1)= alpha/dr^2;
        b(i)= -T(m,i)/dt - qdot/p/Cp;
    end
    A(1,1) = -1; A(1,2) = 1; b(1) = 0;
    A(nx,nx) = 1 + h*dr/k; A(nx,nx-1) = -1; b(nx) = h*dr/k*Tinf;

    T(m+1,:)=inv(A)*b;
end

mesh(r,t,T)
xlabel('x (m)')
ylabel('t (s)')
zlabel('T (^oC)')
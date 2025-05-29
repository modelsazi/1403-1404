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

T = zeros(nt,nx);
T(1,:) = Ti;
T(:,1) = Ti;
T(:,nx)= Ti;

alfa = k/p/Cp;
beta = alfa*dt/dr^2;

for m=1:nt-1
    for i = 2:nx-1
        T(m+1,i) = beta/r(i)*r(i+1)*T(m,i+1)+(1-beta/r(i)*(r(i+1)+r(i)))*T(m,i)+beta*T(m,i-1)+qdot*dt/p/Cp;
    end
    T(m+1,1)  = T(m+1,2);
    T(m+1,nx) = (T(m+1,nx-1) + h*dr*Tinf/k)/(1+h*dr/k); 
end
mesh(r,t,T)
xlabel('r (m)')
ylabel('t (s)')
zlabel('T (^oC)')

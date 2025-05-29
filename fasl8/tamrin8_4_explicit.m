clc;clear;close all

Ri = 0.25; Ro = 0.35; V = 50;

ndx = 100; nx = ndx+1;
dr = (Ro-Ri)/ndx;
r = Ri:dr:Ro;

ndt=5000; Lt = 1000;
nt=ndt+1;
dt=Lt/ndt;
t=0:dt:Lt;

alpha = 10;
Vz = zeros(nt,nx);

for m=1:nt-1
    for i = 2:nx-1
        Vz(m+1,i) = Vz(m,i) + 1/r(i)/alpha*(dr*(Vz(m,i+1)-Vz(m,i))+r(i)*(Vz(m,i+1)-2*Vz(m,i)+Vz(m,i-1)));
    end
    Vz(m+1,1)  = V;
    Vz(m+1,nx) = 0; 
end
mesh(r,t,Vz)
xlabel('r (m)')
ylabel('t (s)')
zlabel('V (m/s)')

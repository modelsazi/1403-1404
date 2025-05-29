clc;clear;close all

Ti = 300; % initial temperature
p = 960; Cp = 2200; R = 0.3; h = 75; qdot = 1000; Tinf = 300;
k0 = 0.16;

ndx = 80; nx = ndx+1;
dr = R/ndx;
r = 0:dr:R;

ndt=100; Lt = 4000;
nt=ndt+1;
dt=Lt/ndt;
t=0:dt:Lt;

alpha = k0*dt/p/Cp/dr^2;
beta = qdot*dt/p/Cp;

T = ones(nt,nx)*Ti;

for m = 2:nt-1
    err = 1;
    T(m+1,:) = T(m,:); % initial guess (last step temperature)
    j = zeros(nx,nx);
    b = zeros(nx,1);
    while err > 1e-3  % for newton-rafson
        for i = 2:nx-1
            j(i,i+1) = -alpha/r(i)*(r(i+1)*(2*T(m+1,i+1)-T(m+1,i)));
            j(i,i) = 1+alpha/r(i) *(r(i+1)*T(m+1,i+1)+r(i)*(2*T(m+1,i)-T(m+1,i-1)));
            j(i,i-1) = -alpha*T(m+1,i);
            b(i) = -(T(m+1,i)-T(m,i)-alpha/r(i)*(r(i+1)*(T(m+1,i+1)^2-T(m+1,i+1)*T(m+1,i))-r(i)*(T(m+1,i)^2-T(m+1,i)*T(m+1,i-1)))-beta);  
        end
        j(1,1) = -1; j(1,2) = 1; b(1) = -(T(m+1,2)-T(m+1,1));
        j(nx,nx) = -k0*2*T(m+1,nx)+k0*T(m+1,nx-1)-h*dr;
        j(nx,nx-1) = k0*T(m+1,nx);
        b(nx) = -(-k0*T(m+1,nx)^2+k0*T(m+1,nx)*T(m+1,nx-1)-h*dr*T(m+1,nx)+h*dr*Tinf);
        dT = j\b;
        T(m+1,:) = T(m+1,:) + dT';
        err = max(abs(dT));
    end
end

mesh(r,t,T)
xlabel('r (m)')
ylabel('t (s)')
zlabel('T (^oC)')
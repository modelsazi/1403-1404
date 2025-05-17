clc; clear; close all

Ea = 15e3; R_gas = 8.314; Ca0 = 1.5; k1 = 200;
q = 0.0063; T = 332.9164; V = 5;

f = @(x,y) q*Ca0/V -q*y/V -k1*exp(-Ea/R_gas/T)*y^2;
x_start = 0;
x_end   = 100;
dx = 0.01;
x = x_start:dx:x_end;
y = zeros(1,length(x));  y(1,1) = 0.0369;
i = 1;
while i < length(x)
    k1 = f(x(i),y(:,i));
    k2 = f(x(i) + dx/2, y(:,i) + k1*dx/2);
    k3 = f(x(i) + dx/2, y(:,i) + k2*dx/2);
    k4 = f(x(i) + dx  , y(:,i) + k3*dx  );
    y(:,i+1)  = y(:,i) + dx/6 * ( k1 + 2*k2 + 2*k3 + k4);
    i = i + 1;
end

plot(x,y(1,:));grid
legend("C_A")
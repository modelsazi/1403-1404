clc; clear; close all

Vx = 20; Ti = 300; density = 800; R = 0.01;
Cp = 420; k = 7e6; L = 50; 

f = @(x,y) (2*k)/(Vx*density*Cp*R*y);
x_start = 0;
x_end   = L;
dx = 0.01;
x = x_start:dx:x_end;
y = zeros(1,length(x));  y(1,1) = Ti;
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
legend("Temp")
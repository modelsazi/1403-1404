clc; clear; close all

h = 85; R = 0.2; rho = 1150; Tb = 400;
f = @(x,y) 2*h/(R*rho)*(Tb-y)/y^0.5;
x_start = 0;
x_end   = 250;
dx = 0.001;
x = x_start:dx:x_end;
y = zeros(1,length(x));  y(1) = 298;
i = 1;
while i < length(x)
    k1 = f(x(i),y(:,i));
    k2 = f(x(i) + dx/2, y(:,i) + k1*dx/2);
    k3 = f(x(i) + dx/2, y(:,i) + k2*dx/2);
    k4 = f(x(i) + dx  , y(:,i) + k3*dx  );
    y(:,i+1)  = y(:,i) + dx/6 * ( k1 + 2*k2 + 2*k3 + k4);
    i = i + 1;
end
%figure(2)
plot(x,y(1,:));grid
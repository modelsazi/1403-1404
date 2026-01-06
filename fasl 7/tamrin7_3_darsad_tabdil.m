clc; clear; close all

k1 = 0.03; k2 = 0.01; B0 = 1.7; A0 = 1;
f = @(x,y) k1*(1-y(1))*(B0-A0*y(1)) - k2*y(1)^2*A0;

x_start = 0;
x_end   = 220;
dx = 0.01;
x = x_start:dx:x_end;
y = zeros(1,length(x));  y(1,1) = 0;

%% runge-kutta fourth order method
i = 1;
while i < length(x)
    k1 = f(x(i),y(:,i));
    k2 = f(x(i) + dx/2, y(:,i) + k1*dx/2);
    k3 = f(x(i) + dx/2, y(:,i) + k2*dx/2);
    k4 = f(x(i) + dx  , y(:,i) + k3*dx  );
    y(:,i+1)  = y(:,i) + dx/6 * ( k1 + 2*k2 + 2*k3 + k4);
    i = i + 1;
end
plot(x,y(1,:)); grid
legend(["p"])
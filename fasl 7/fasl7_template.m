clc; clear; close all

%% defining constants and equations (changes for every question)
f = @(x,y) sin(x)/y + y;
x_start = 0;
x_end   = 4;
dx = 0.001;
x = x_start:dx:x_end;
y = zeros(1,length(x));  y(1) = 1;

%% Euler method
i = 1;
while i < length(x)
    k1 = f(x(i),y(:,i));
    y(:,i+1) = y(:,i) + k1*dx;
    i = i + 1;
end
figure(1)
plot(x,y(1,:));grid
title("Euler method")

%% Heon method
%i = 1;
%while i < length(x)
%    k1 = f(x(i),y(:,i));
%    k2 = f(x(i) + dx, y(:,i) + k1*dx);
%    y(:,i+1)  = y(:,i) + dx/2 * ( k1 + k2);
%    i = i + 1;
%end
%figure(2)
%plot(x,y(1,:));grid
%title("heon method")

%% runge-kutta third order method
%i = 1;
%while i < length(x)
%    k1 = f(x(i),y(:,i));
%    k2 = f(x(i) + dx/2, y(:,i) + k1*dx/2);
%    k3 = f(x(i) + dx, y(:,i) + 2*k2*dx - dx*k1);
%    y(:,i+1)  = y(:,i) + dx/6 * ( k1 + 4*k2 + k3);
%    i = i + 1;
%end
%figure(3)
%plot(x,y(1,:));grid
%title("runge-kutta third order")

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
figure(4)
plot(x,y(1,:));grid
title("runge-kutta fourth order")

%% runge-kutta merson method
%i = 1;
%while i < length(x)
%    k1 = f(x(i),y(:,i));
%    k2 = f(x(i) + dx/3, y(:,i) + k1*dx/3);
%    k3 = f(x(i) + dx/3, y(:,i) + k1*dx/6 + k2*dx/6);
%    k4 = f(x(i) + dx/2, y(:,i) + k1*dx/8 + 3*k3*dx/8);
%    k5 = f(x(i) + dx  , y(:,i) + k1*dx/2 - 3*k3*dx/2 + 2*dx*k4);
%    y(:,i+1)  = y(:,i) + dx/6 * ( k1 + 4*k4 + k5);
%    i = i + 1;
%end
%figure(5)
%plot(x,y(1,:));grid
%title("runge-kutta merson")

%% Gill method
i = 1;
while i < length(x)
    k1 = f(x(i),y(:,i));
    k2 = f(x(i) + dx/2, y(:,i) + k1*dx/2);
    k3 = f(x(i) + dx/2, y(:,i) + (-1/2 + 1/sqrt(2))*dx*k1 + (1/2 - 1/sqrt(2))*dx*k2);
    k4 = f(x(i) + dx  , y(:,i) - k1*dx/sqrt(2) + (1/2 + 1/sqrt(2))*dx*k2);
    y(:,i+1)  = y(:,i) + dx/6 * ( k1 + 2*(1-1/sqrt(2))*k2 + 2*(1+1/sqrt(2))*k3 + k4);
    i = i + 1;
end
figure(6)
plot(x,y(1,:));grid
title("Gill method")

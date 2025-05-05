clc; clear; close all

f = @(x,y) sin(x)/y + y;

dx = 0.1;
X = 0:dx:4  ; % 0,0.1,0.2, .... , 4  
Y = ones(1,length(X));
m = 1;

while X(m) < 4
    k1 = f(X(m),Y(m));
    k2 = f(X(m) + dx, Y(m) + k1*dx);
    Y(m+1)  = Y(m) + dx/2 * ( k1 + k2);
    m = m + 1; 
end

plot(X,Y);grid

clc; clear; close all

k = 0.018; L = 0.05; h = 15; Tinf=298;

dx = 0.001; x = 0:dx:L; n = length(x); % solution A
A = zeros(n); B = zeros(n,1);

for i = 2:n - 1
    A(i,i-1) = x(i)^4;
    A(i,i) = -(x(i+1)^4 + x(i)^4 + dx^2*x(i)^2/2/k);
    A(i,i+1) = x(i+1)^4;
    B(i) = -dx^2/2/k*x(i)^2*Tinf;
end

A(1,1) = 1; B(1) = Tinf;
A(n,n)= 1 ; B(n) = 373;
T = A\B;
plot(x,T); grid

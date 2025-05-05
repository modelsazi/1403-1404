clc; clear; close all

qdot = 1000; k = 0.16; R = 0.3; h = 85; Tinf = 300;
dr = 0.01; r = 0:dr:R; n = length(r); % solution A
%n = 100;    dr = (R-0)/n;    r = @(i) (i-1)*dr;  % solution B
A = zeros(n); B = zeros(n,1);

for i = 2:n - 1
    A(i,i-1) = r(i);
    A(i,i) = -(r(i+1) + r(i));
    A(i,i+1) = r(i+1);
    B(i) = -r(i)* dr^2*qdot/k;
end

A(1,1) = -1; A(1,2) = 1;
A(n,n)= (k+h*dr); A(n,n-1) = -k; B(n) = h*dr*Tinf;
T = A\B;
% solution A
plot(r,T); grid
% solution B
%plot(r(1:n),T); grid
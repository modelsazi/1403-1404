clc; clear; close all

%% Alef
L = 1; D = 1; u = 5; k = 10;  
nx = 100; n = nx + 1;  dz = (L-0)/nx; z = 0:dz:L; 
C0 = 5;
j = zeros(n); b = zeros(n,1);
C = ones(n,1);  err = 1;

while err > 1e-4 % newton-rafson loop
    for i = 2:n - 1  % calculating coefficients
        j(i,i-1) = D;
        j(i,i)   = -2*D+u*dz-2*k*dz^2*C(i);
        j(i,i+1) = D-u*dz;
        b(i)     = -((D-u*dz)*C(i+1)+(-2*D+u*dz)*C(i)-k*dz^2*C(i)^2+D*C(i-1));
    end
    j(1,1) = -(D+u*dz); j(1,2) = D; b(1) = -(D*C(2)-(D+u*dz)*C(1)+u*dz*C0);
    j(n,n)= 1; j(n,n-1) = -1; b(n) = -(C(n) - C(n-1));
    dC = j\b;
    err = max(abs(dC));
    C = C + dC;
end
subplot(1,2,1)
plot(z,C); grid
title("alef")

%% B
L = 1; D = 1; u = 5; k = 10;  
nx = 100; n = nx + 1;  dz = (L-0)/nx; z = 0:dz:L; 
C0 = 5;
A = zeros(n); B = zeros(n,1);

for i = 2:n-1
    A(i,i-1) = D;
    A(i,i) = -(2*D+u*dz-k*dz^2);
    A(i,i+1) = D-u*dz;
    B(i) = 0;
end

A(1,1)= -(D+u*dz); A(1,2) = D; B(1) = -u*dz*C0;
A(n,n)= 1 ; A(n,n-1) = -1; B(n) = 0;
C = A\B;

subplot(1,2,2)
plot(z,C); grid
title("B")
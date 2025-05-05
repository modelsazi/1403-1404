clc; clear; close all

k = 0.16; R = 0.3; h = 85; Tinf = 300;
dr = 0.01; r = 0:dr:R; n = length(r); % solution A
%n = 100;    dr = (R-0)/n;    r = @(i) (i-1)*dr;  % solution B
j = zeros(n); b = zeros(n,1);
T = ones(n,1) * 350;  err = 1;
%m = 1;
while err > 1e-4 % newton-rafson loop
    for i = 2:n - 1  % calculating coefficients
        j(i,i-1) = 0.16*T(i)*r(i);
        j(i,i)   = -0.16*r(i+1)*T(i+1)-0.16*2*r(i)*T(i)+0.16*r(i)*T(i-1)+100*r(i)*dr^2*1.4*(T(i)-Tinf)^0.4;
        j(i,i+1) = 0.16*2*r(i+1)*T(i+1)-0.16*r(i+1)*T(i);
        b(i)     = -(0.16*r(i+1)*T(i+1)^2-0.16*r(i+1)*T(i+1)*T(i)-0.16*r(i)*T(i)^2+0.16*r(i)*T(i)*T(i-1)+100*r(i)*dr^2*(T(i)-Tinf)^1.4);
    end
    j(1,1) = -1; j(1,2) = 1; b(1) = -(T(2) - T(1));
    j(n,n)= 0.16*2*T(n)+h*dr-0.16*T(n-1); j(n,n-1) = -0.16*T(n); b(n) = -(0.16*T(n)^2+h*dr*T(n)-0.16*T(n)*T(n-1)- h*dr*Tinf);
    dx = j\b;
    err = max(abs(dx));
    T = T + dx;
    %m = m + 1;
end
% solution A
plot(r,T); grid
%disp(m)
% solution B
%plot(r(1:n),T); grid





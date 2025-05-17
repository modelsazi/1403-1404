clc; clear; close all

%% Alef
L = 1; landa = 1; b = 4;
dx = 0.001; x = 0:dx:L; n = length(x);
A = zeros(n); B = zeros(n,1);

for i = 2:n - 1
    A(i,i-1) = x(i)^2;
    A(i,i) = -(2*x(i)*dx+2*x(i)^2+landa^2*dx^2*x(i));
    A(i,i+1) = (2*x(i)*dx+x(i)^2);
    B(i) = -landa^2*x(i)*dx^2;
end

A(1,1) = -1; A(1,2) = 1; B(1) = 0;
A(n,n)= 1 ; B(n) = b;
y = A\B;
subplot(2,2,1)
plot(x,y); grid
title("Alef")
%% B
f = @(x,y) [y(2)
            -2/x*y(2)+landa^2*y(1)/x+landa^2/x];
x_start = 0.1;
x_end   = L;
dx = 0.001;
x = x_start:dx:x_end;
y = zeros(2,length(x));  y(1,1) = b; y(2,1) = 0;
i = 1;
while i < length(x)
    k1 = f(x(i),y(:,i));
    k2 = f(x(i) + dx/2, y(:,i) + k1*dx/2);
    k3 = f(x(i) + dx/2, y(:,i) + k2*dx/2);
    k4 = f(x(i) + dx  , y(:,i) + k3*dx  );
    y(:,i+1)  = y(:,i) + dx/6 * ( k1 + 2*k2 + 2*k3 + k4);
    i = i + 1;
end
subplot(2,2,2)
plot(x,y(1,:),x,y(2,:));grid
legend("y","p")
title("B")

%% C
f = @(x,y) [y(2)
            -2/y(1)*y(2)+landa^2*x/y(1)+landa^2*x/y(1)^2];
x_start = 0.1;
x_end   = L;
dx = 0.001;
x = x_start:dx:x_end;
y = zeros(2,length(x));  y(1,1) = b; y(2,1) = 0;
i = 1;
while i < length(x)
    k1 = f(x(i),y(:,i));
    k2 = f(x(i) + dx/2, y(:,i) + k1*dx/2);
    k3 = f(x(i) + dx/2, y(:,i) + k2*dx/2);
    k4 = f(x(i) + dx  , y(:,i) + k3*dx  );
    y(:,i+1)  = y(:,i) + dx/6 * ( k1 + 2*k2 + 2*k3 + k4);
    i = i + 1;
end
subplot(2,2,3)
plot(x,y(1,:),x,y(2,:));grid
legend("y","p")
title("C")

%% D
L = 1; landa = 1; b = 4;
dx = 0.001; x = 0:dx:L; n = length(x); 
j = zeros(n,n); B = zeros(n,1);
y = ones(1,length(x)); err = 1;
while err > 1e-4 % newton-rafson loop
    for i = 2:n - 1  % calculating coefficients
        j(i,i-1) = y(i)^2;
        j(i,i)   = 2*y(i+1)*dx-4*y(i)*dx+2*y(i)*y(i+1)-6*y(i)^2+2*y(i)*y(i-1)-landa^2*x(i)*dx^2;
        j(i,i+1) = 2*y(i)*dx+y(i)^2;
        B(i)     = -(2*y(i)*y(i+1)*dx-2*y(i)^2*dx+y(i)^2*y(i+1)-2*y(i)^3+y(i)^2*y(i-1)-landa^2*x(i)*y(i)*dx^2+landa^2*x(i)*dx^2);
    end
    j(1,1) = -1; j(1,2) = 1; B(1) = -(y(2) - y(1));
    j(n,n)= 1; B(n) = b-y(n);
    dy = j\B;
    err = max(abs(dy));
    y = y + dy;
end
subplot(2,2,4)
plot(x,y); grid
title("D")
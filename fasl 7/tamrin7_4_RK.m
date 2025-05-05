clc; clear; close all

f = @(x,y) [3+2/3*y(1)-5/3*y(1)-3*y(1)*y(2)*1+2*1*y(3)^2
           2+2/3*y(2)-5/3*y(2)-3*y(1)*y(2)*1+2*1*y(3)^2
           2/3*y(3)-5/3*y(3)+3*y(1)*y(2)*1-2*1*y(3)^2];
x_start = 0;
x_end   = 40;
dx = 0.001;
x = x_start:dx:x_end;
y = zeros(3,length(x));  y(1) = 1; y(2) = 1; y(3) = 1;
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
plot(x,y(1,:),x,y(2,:),x,y(3,:));grid
legend("A","B","C")
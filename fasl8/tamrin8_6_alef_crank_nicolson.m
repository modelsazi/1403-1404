clc;clear;close all

ui = 1; L = 1;

nx = 40; ndx = nx - 1;
dx = L/ndx;
x = 0:dx:L;

nbeta= 100; ndbeta = nbeta - 1;
Lbeta = 2;
dbeta=Lbeta/ndbeta;
beta=0:dbeta:Lbeta;

u=zeros(nbeta,nx);
u(1,:)=ui;

for m = 1:nbeta-1
    A = zeros(nx,nx); b = zeros(nx,1);
    for i = 2:nx-1
        A(i,i-1)= 1/(4*dx^2);
        A(i,i)=  -(1/(2*dx^2)+1/dbeta);
        A(i,i+1)= 1/(4*dx^2);
        b(i)= -u(m,i)/dbeta-1/4*(u(m,i+1)-2*u(m,i)+u(m,i-1))/(dx^2);
    end
    A(1,1) = 1; b(1) = 0;
    A(nx,nx) = 1; A(nx,nx-1) = -1; b(nx) = 0;
    u(m+1,:) = (A\b);
end
mesh(x,beta,u)
xlabel('x (m)')
ylabel('\beta (s)')
zlabel('u (^oC)')
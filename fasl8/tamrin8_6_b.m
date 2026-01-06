clc;clear;close all

nx = 30;
ndx=nx-1;
dx=1/ndx;
x=0:dx:1;

ny = 30;
ndy=ny-1;
dy=1/ndy;
y=0:dy:1;

n=nx*ny;

A=zeros(n,n);
b=zeros(n,1);

for k=ny+2:n-ny
    A(k,k+ny)=1/dx^2;
    A(k,k-ny)=1/dx^2;
    A(k,k+1)=1/(dy^2);
    A(k,k-1)=1/(dy^2);
    A(k,k)= -(2/(dx^2)+2/(dy^2));
    b(k)=-2;
end

for k=1:ny  % left
    A(k,:)=0; A(k,k)= 1; b(k)= 0;
end

for k=1:ny:n-ndy  % bottom
    A(k,:)=0; A(k,k+1)=1; A(k,k)=-1; b(k)=0;
end

for k=n-ndy:n % right
    A(k,:)=0; A(k,k)= 1; b(k)= 0;
end

for k=ny:ny:n  % up
    A(k,:)=0; A(k,k)=1; A(k,k-1)=-1; b(k)=dy;
end

u=A\b;
u1 = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        u1(i,j)=u((i-1)*ny+j);
    end
end

surf(y,x,u1)
xlabel('y (m)')
ylabel('x (m)')
zlabel('u')
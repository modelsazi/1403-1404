clc;clear;close all

w=0.1; H=0.1; kk=15; Ta=100; qb=900; Tinf=25; h=30; qdot=10000;

ndx=30;
nx=ndx+1;
dx=w/ndx;
x=0:dx:w;

ndy=30;
ny=ndy+1;
dy=H/ndy;
y=0:dy:H;

n=nx*ny;

A=zeros(n,n);
b=zeros(n,1);

for k=ny+1:n-ny
    A(k,k+ny)=1;
    A(k,k-ny)=1;
    A(k,k+1)=1;
    A(k,k-1)=1;
    A(k,k)=-4;
    b(k)=-qdot*dx^2/kk;
end

for k=1:ny  % left 
    A(k,:)=0; A(k,k+ny)=1; A(k,k)=-1; b(k)=0;
end

for k=1:ny:n-ndy  % bottom
    A(k,:)=0; A(k,k+1)=1; A(k,k)=-1; b(k)=-qb*dy/kk;
end

for k=n-ndy:n % right
    A(k,:)=0; A(k,k-ny)=-1; A(k,k)=1+h*dx/kk; b(k)=h*dx*Tinf/kk;
end

for k=ny:ny:n  % up
    A(k,:)=0; A(k,k)=1; b(k)=Ta*(1+sin(pi*(k-ny)/ny*dx/w)/4); 
end

T=A\b;
T1 = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        T1(i,j)=T((i-1)*ny+j);
    end
end

surf(y,x,T1)
xlabel('y (m)')
ylabel('x (m)')
zlabel('T (^oC)')

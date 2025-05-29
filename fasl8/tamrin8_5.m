clc;clear;close all

kk=0.34; Ta1=100; Ta2 = 20; Tinf = 120; h=15; 
Ri = 0.05; Ro = 0.35; THETA = pi;

ndr=30;
nr=ndr+1;
dr=(Ro-Ri)/ndr;
r=Ri:dr:Ro;

ndy=30; % y = theta
ny=ndy+1;
dy=THETA/ndy;
y=0:dy:THETA;

n=nr*ny;

A=zeros(n,n);
b=zeros(n,1);

for k=ny+1:n-ny
    A(k,k+ny)=r(fix(k/ny)+1)/dr^2+1/2/dr;
    A(k,k-ny)=(-1/2/dr+r(fix(k/ny)+1)/dr^2);
    A(k,k+1)=1/dy^2;
    A(k,k-1)=1/dy^2;
    A(k,k)=(-2/dy^2-2*r(fix(k/ny)+1)/dr^2);
    b(k)=0;
end

for k=1:ny  % left 
    A(k,:)=0; A(k,k+ny)= 1; A(k,k)=-1; b(k)= 0;
end

for k=1:ny:n-ndy  % bottom
    A(k,:)=0; A(k,k)=1; b(k)=Ta1;
end

for k=n-ndy:n % right
    A(k,:)=0; A(k,k-ny)=+kk/dr; A(k,k)=(-kk/dr-h); b(k)=-h*Tinf;
end

for k=ny:ny:n  % up
    A(k,:)=0; A(k,k)=1; b(k)=Ta2; 
end

T=A\b;
T1 = zeros(nr,ny);
for i=1:nr
    for j=1:ny
        T1(i,j)=T((i-1)*ny+j);
    end
end

surf(y,r,T1)
xlabel('y (m)')
ylabel('r (m)')
zlabel('T (^oC)')
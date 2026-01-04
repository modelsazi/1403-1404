clc;clear;close all

kk=0.34; Ta1=20; Ta2 = 100; Tinf = 30; h=15; 
Ri = 0.05; Ro = 0.35; THETA = pi;

ndr=30;
nr=ndr+1;
dr=(Ro-Ri)/ndr;
r=Ri:dr:Ro; 

ndt=30; % y = theta
nt=ndt+1;
dt=THETA/ndt;
t=0:dt:THETA;

n=nr*nt;

A=zeros(n,n);
b=zeros(n,1);

for k=nt+1:n-nt
    rk = r(fix(k/nt)+1);
    A(k,k+nt)=rk/dr^2+1/(2*dr);
    A(k,k-nt)=(-1/2/dr+rk/dr^2);
    A(k,k+1)=1/rk/dt^2;
    A(k,k-1)=1/rk/dt^2;
    A(k,k)=(-2/rk/dt^2-2*rk/dr^2);
    b(k)=0;
end

for k=1:nt  % left
    A(k,:)=0; A(k,k+nt)= 1; A(k,k)=-1; b(k)= 0;
end

for k=1:nt:n-ndt  % bottom
    A(k,:)=0; A(k,k)=1; b(k)=Ta2;
end

for k=n-ndt:n % right
    A(k,:)=0; A(k,k-nt)=+kk/dr; A(k,k)=(-kk/dr-h); b(k)=-h*Tinf;
end

for k=nt:nt:n  % up
    A(k,:)=0; A(k,k)=1; b(k)=Ta1; 
end

T=A\b;
T1 = zeros(nr,nt);
for i=1:nr
    for j=1:nt
        T1(i,j)=T((i-1)*nt+j);
    end
end

%[r,t] = meshgrid(t,r);
%[x,y] = pol2cart(r,t);
%surf(y,x,T1)
surf(t,r,T1)
xlabel('theta (m)')
ylabel('r (m)')
zlabel('T (^oC)')
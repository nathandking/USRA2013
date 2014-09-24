function [uu]=intdiff(u,ht,hx,n,aij)

% compute summation for uhat.
summ=zeros(n,1);
summ(2:n-1)=aij(2:n-1,2:n)*u(2:n)+aij(2:n-1,1:n-1)*u(1:n-1);

uhat=zeros(n,1);
uhat(2:n-1)=u(2:n-1)-0.25*ht*hx*summ(2:n-1);

% compute summation for u.
summm=zeros(n,1);
summm(2:n-1)=aij(2:n-1,2:n)*uhat(2:n)+aij(2:n-1,1:n-1)*uhat(1:n-1);

uu=zeros(n,1);
uu(2:n-1)=u(2:n-1)-0.5*ht*hx*summm(2:n-1);

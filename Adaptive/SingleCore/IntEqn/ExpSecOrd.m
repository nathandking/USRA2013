function [yhat]=ExpSecOrd(y,h,n,hx)

k1=zeros(n,1);
k1(2:n-1)=y(2:n-1)+0.5*h*(y(1:n-2)-2*y(2:n-1)+y(3:n))/(hx^2);

yhat=zeros(n,1);
yhat(2:n-1)=y(2:n-1)+h*(k1(1:n-2)-2*k1(2:n-1)+k1(3:n))/(hx^2);
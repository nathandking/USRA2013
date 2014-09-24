function [yhat]=ExpSecOrd(y,h,n,hx)

y1=zeros(n,1);
y1(2:n-1)=y(2:n-1)+0.5*h*(y(1:n-2)-2*y(2:n-1)+y(3:n))/(hx^2);

yhat=zeros(n,1);
yhat(2:n-1)=y(2:n-1)+h*(y1(1:n-2)-2*y1(2:n-1)+y1(3:n))/(hx^2);

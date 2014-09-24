function [yy]=ExpMid(y,h,hx,n)

y1=zeros(n,1);
y1(2:n-1)=y(2:n-1)+(0.5*h/hx^2)*(y(1:n-2)-2*y(2:n-1)+y(3:n));

yy=zeros(n,1);
yy(2:n-1)=y(2:n-1)+(h/hx^2)*(y1(1:n-2)-2*y1(2:n-1)+y1(3:n));
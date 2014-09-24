function [xx]=g(x,y,h,n,hx)

xx(1)=x(1);
xx(2:n-1)=x(2:n-1)-y(2:n-1)-h*(x(1:n-2)-...
    2*x(2:n-1)+x(3:n))/(2*hx^2);
xx(n)=x(n);

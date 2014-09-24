function [gz]=trap(x,y,h,n,hx,delta,aij)

gz=x-(y+0.5*h*(f(h,y,n,hx,delta,aij)+f(h,x,n,hx,delta,aij)));
gz(1)=x(1);
gz(n)=x(n);
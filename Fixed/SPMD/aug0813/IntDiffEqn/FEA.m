function [xx]=FEA(x,y,h,hx,n,delta,aij)

xx=x-y-h*f(h,x,n,hx,delta,aij);
xx(1)=x(1);
xx(n)=x(n);
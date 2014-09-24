function [gz]=trap(x,y,h,n,hx,delta,alpha,p)
%-----------------------------------------------------------------------------
%    The standard trapezoidal rule for differential equations. This function
%    is written to be used by fsolve to solve the nonlinear implicit 
%    equations.
%-----------------------------------------------------------------------------

gz=x-(y+0.5*h*(fbeck(h,y,n,hx,delta,alpha,p)+fbeck(h,x,n,hx,delta,alpha,p)));
gz(1)=x(1);
gz(n)=x(n);

function [psi]=PsiFEA(x,y,h,n,delta,hx)
%--------------------------------------------------------------------------
%   This function defines the implicit function to be solved by fsolve.
%   This function corresponds to the first order backward Euler B-method
%   with F(u)=exp(u).
%--------------------------------------------------------------------------

psi(1)=x(1);
psi(2:n-1)=x(2:n-1)-(y(2:n-1)+h*(x(1:n-2)-2*x(2:n-1)+x(3:n))/hx^2+delta*h*exp(x(2:n-1)));
psi(n)=x(n);
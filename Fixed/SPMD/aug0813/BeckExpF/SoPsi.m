function [psi]=SoPsi(x,y,h,n,delta,hx)
%--------------------------------------------------------------------------
%   This function defines the implicit function to be solved by fsolve.
%   This function corresponds to the second order forward Euler B-method
%   with F(u)=exp(u).
%--------------------------------------------------------------------------

psi(1)=x(1);
psi(2:n-1)=x(2:n-1)-(0.5*h/hx^2)*(x(1:n-2)-2*x(2:n-1)+x(3:n))+...
    log(exp(-y(2:n-1))-0.5*delta*h);
psi(n)=x(n);
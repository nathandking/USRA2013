function [psi]=SoPsi(x,y,h)
%--------------------------------------------------------------------------
%   This function defines the implicit function to be solved by fsolve.
%   This function corresponds to the second order forward Euler B-method
%   with F(u)=exp(u).
%--------------------------------------------------------------------------

global n hx;

psi(1)=x(1);
psi(2:n-1)=x(2:n-1)-y(2:n-1)-(h/(hx^2))*(x(1:n-2)-2*x(2:n-1)+x(3:n));
psi(n)=x(n);


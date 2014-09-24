function [psi]=BE(x,y,k)
%--------------------------------------------------------------------------
%   This function defines the implicit function to be solved by fsolve.
%   This function corresponds to the first order backward Euler method.
%--------------------------------------------------------------------------

global n hx alpha;

psi(1)=x(1);
psi(2:n-1)=x(2:n-1)-(y(2:n-1)+k*(((x(1:n-2)-2*x(2:n-1)+x(3:n))/(hx^2))+...
    (x(2:n-1)).^(1+alpha)));
psi(n)=x(n);
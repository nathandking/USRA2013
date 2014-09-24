function [psi]=SoPsi(x,y,h,n,hx,delta,alpha,p)
%--------------------------------------------------------------------------
%   This function defines the implicit function to be solved by fsolve.
%   This function corresponds to the second order forward Euler B-method
%   with F(u)=(u+alpha)^p.
%--------------------------------------------------------------------------

psi(1)=x(1);
psi(2:n-1)=x(2:n-1)-(0.5*h/hx^2)*(x(1:n-2)-2*x(2:n-1)+x(3:n))-...
    ((y(2:n-1)+alpha).^(1-p)+0.5*(1-p)*delta*h).^(1/(1-p))+alpha;
psi(n)=x(n);
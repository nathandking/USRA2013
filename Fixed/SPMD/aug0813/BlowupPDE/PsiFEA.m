function [psi]=PsiFEA(x,y,h,n,hx,delta,alpha,p)
%--------------------------------------------------------------------------
%   This function defines the implicit function to be solved by fsolve.
%   This function corresponds to the first order backward Euler method
%   with F(u)=(u+alpha)^p.
%--------------------------------------------------------------------------

psi(1)=x(1);
psi(2:n-1)=x(2:n-1)-(y(2:n-1)+h*(((x(1:n-2)-2*x(2:n-1)+x(3:n))/(hx^2))+...
    delta*(x(2:n-1)+alpha*ones(n-2,1)).^p));
psi(n)=x(n);
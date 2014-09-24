function [psi]=CN(x,y,k,h,alpha)
%--------------------------------------------------------------------------
%   This function defines the implicit function to be solved by fsolve.
%   This function corresponds to the second order Crank-Nicolson method for
%   u_t=u_xx+u^(1+\alpha).
%--------------------------------------------------------------------------
n=length(y);

psi(1)=x(1);
psi(2:n-1)=x(2:n-1)-(y(2:n-1)+0.5*k*(((x(1:n-2)-2*x(2:n-1)+x(3:n))/(h^2))+...
    (x(2:n-1)).^(1+alpha)+((y(1:n-2)-2*y(2:n-1)+y(3:n))/(h^2))+...
    (y(2:n-1)).^(1+alpha)));
psi(n)=x(n);
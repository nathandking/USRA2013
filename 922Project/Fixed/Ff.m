function [dy]=f(~,y)
%--------------------------------------------------------------------------
%   This function is used by ode15s to solve a ODE system obtained by 
%   applying the MOL to 
%                   $$u_t=u_{xx}+F(u),$$
%   with Dirichlet boundary conditions and $F(u)=u^{1+\alpha}$.
%--------------------------------------------------------------------------
global n hx alpha;

dy=zeros(n,1);
dy(2:n-1)=(y(1:n-2)-2*y(2:n-1)+y(3:n))/(hx^2)+y(2:n-1).^(1+alpha);

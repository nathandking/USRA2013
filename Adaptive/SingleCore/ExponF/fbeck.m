function [dy]=fbeck(~,y,n,hx,delta)
%--------------------------------------------------------------------------
%   This function is used by ode45 to solve a ODE system obtained by 
%   applying the MOL to 
%                   u_t=u_xx+delta*F(u),
%   with Dirichlet boundary conditions and F(u)=exp(u).
%--------------------------------------------------------------------------

dy=zeros(n,1);
dy(2:n-1)=((y(1:n-2)-2*y(2:n-1)+y(3:n))/(hx^2))+delta*exp(y(2:n-1));



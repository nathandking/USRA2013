function [gz]=trap(x,y,h,lambda,p)
%--------------------------------------------------------------------------
%   This function is used be fzero to solve the implicit equation of the
%   trapezoidal rule applied to y'=lambda*y+y^p.
%--------------------------------------------------------------------------

gz=x-(y+0.5*h*(lambda*y+y^p+lambda*x+x^p));
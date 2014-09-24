function [gx]=g(x,y,h,p,lambda)
%--------------------------------------------------------------------------
%   Function used for fzero to solve the implicit backward Euler method.
%--------------------------------------------------------------------------

gx=x-(y+h*(lambda*x+x^p));
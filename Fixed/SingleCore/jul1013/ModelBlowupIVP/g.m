function [gx]=g(x,y,h)
global lambda p;

gx=x-(y+h*(lambda*x+x^p));
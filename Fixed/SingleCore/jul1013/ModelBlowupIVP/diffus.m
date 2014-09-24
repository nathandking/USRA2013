function [y]=diffus(y0,h)
global lambda;

k1=y0+h*lambda*y0;
y=y0+(h/2)*(lambda*y0+lambda*k1);

%y=y0+h*lambda*y0;
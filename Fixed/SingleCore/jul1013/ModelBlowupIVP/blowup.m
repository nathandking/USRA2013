function [y]=blowup(y0,t)
global p ;

y=((1-p)*t+y0.^(1-p)).^(1./(1-p));
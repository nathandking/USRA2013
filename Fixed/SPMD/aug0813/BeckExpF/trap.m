function [gz]=trap(x,y,h,n,hx,delta)

gz=x-(y+0.5*h*(fbeck(h,y,n,hx,delta)+fbeck(h,x,n,hx,delta)));
gz(1)=x(1);
gz(n)=x(n);
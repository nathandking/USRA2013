function [yhat]=ImpSecOrd(y,h,n,hx)

options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,'Display','off');
k1=fsolve(@(x) g(x,y,h,n,hx),y,options);

yhat=zeros(n,1);
yhat(2:n-1)=y(2:n-1)+0.5*h*(k1(1:n-2)-2*k1(2:n-1)+k1(3:n))/(hx^2);
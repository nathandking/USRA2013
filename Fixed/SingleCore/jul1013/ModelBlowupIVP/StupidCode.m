function [ge, geFEA]=StupidCode(ht,t0,tf,y0)
global p lambda lamp;

z(1)=y0;
y(1)=y0;
t=t0:ht:tf;
    
   for i=1:length(t)
    % Integrate one step.
    z(i+1)=fzero(@(x) g(x,z(i),ht),z(i)); 
    
    k=1;
    x(1)=y(i);
    myerr=1;
    while myerr~=0 && k<100
    x(k+1)=x(k)-((((1-ht*lambda)*x(k))-(y(i)^(1-p)+(1-p)*ht)^(1/(1-p)))/(1-ht*lambda)); 
      if norm(x(k+1)-x(k),2) < 1e-12
                myerr =0;
            end
            k=k+1;
    end
    y(i+1)=x(k);
    
    yexact(i)=(((y0^(1-p)+(1/lambda))*exp(lamp*t(i)))-(1/lambda))^(1/(1-p));
    end
ge=abs(y(i)-yexact(i)');
geFEA=abs(z(i-1)-yexact(i)');
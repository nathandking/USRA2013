function GE=BlowupODEBackwardEuler(t0,tf,y0,h,p,lambda,yexact)
%--------------------------------------------------------------------------
%   This function carries out the integration of the IVP
%               y'=lambda*y+y^p,  y(0)=y0,                            (1)
%   with FEA and SpFEA.
%--------------------------------------------------------------------------

   ht=h(labindex);      % unique step size corresponding to worker index. 
   t=t0:ht:tf;
   y=zeros(length(t),1);
   z=zeros(length(t),1);
   y(1)=y0;
   z(1)=y0;
    
   for i=1:length(t);
       % Integrate one step.
       z(i+1)=fzero(@(x) g(x,z(i),ht,p,lambda),z(i));    % FEA.
       
       k=1;                 % SpFEA.
       x=zeros(50,1);
       x(1)=y(i);
       myerr=1;
       while myerr~=0 && k<50;
           x(k+1)=x(k)-((((1-ht*lambda)*x(k))-(y(i)^(1-p)+(1-p)*ht)^...
               (1/(1-p)))/(1-ht*lambda)); 
           if norm(x(k+1)-x(k),2) < 1e-12
               myerr =0;
           end
           k=k+1;
       end
       y(i+1)=x(k);
   end
   % calculate absolute errors.
   geA=abs(y(i)-yexact);
   geFEA=abs(z(i)-yexact);
   GE=[geA,geFEA];



function GE=SoMethods(t0,tf,y0,h,p,lambda,yexact)
%--------------------------------------------------------------------------
%   This function carries out the integration of the IVP
%               y'=lambda*y+y^p,  y(0)=y0,                            (1)
%   with FE and SpFE.
%--------------------------------------------------------------------------

    ht=h(labindex);      % unique step size corresponding to worker index. 
    t=t0:ht:tf;
    y=zeros(length(t),1);
    ybeck=zeros(length(t),1);
    z=zeros(length(t),1);
    w=zeros(length(t),1);
    y(1)=y0;
    ybeck(1)=y0;
    z(1)=y0;
    w(1)=y0;
    
    for i=1:length(t)
        % integrate one step for all four methods.
        z1=z(i)+0.5*ht*(lambda*z(i)+z(i)^p);                  % EMR.
        z(i+1)=z(i)+ht*(lambda*z1+z1^p);
        
        w(i+1)=fzero(@(x) trap(x,w(i),ht,lambda,p),w(i));     % TR.
    
        ybeck(i+1)=((1-p)*0.5*ht+(((1+0.5*lambda*ht)/(1-...   % SoSpFE
            .5*lambda*ht))^(1-p))*(0.5*(1-p)*ht+ybeck(i)^(1-p)))^(1/(1-p));
        
        y1=diffus(y(i),ht/2,lambda);                                 %Strang.
        y2=blowup(y1,ht,p);
        y(i+1)=diffus(y2,ht/2,lambda);
    end
    % calculate absolute errors.
    ge=abs(ybeck(i)-yexact);
    geStrang=abs(y(i)-yexact);
    geExpMid=abs(z(i)-yexact);
    geTrap=abs(w(i)-yexact);
    GE=[ge, geStrang, geExpMid, geTrap];
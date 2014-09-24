function GE=BlowupODEForwardEuler(t0,tf,y0,h,p,lambda,yexact)
%--------------------------------------------------------------------------
%   This function carries out the integration of the IVP
%               y'=lambda*y+y^p,  y(0)=y0,                            (1)
%   with FE and SpFE.
%--------------------------------------------------------------------------

    ht=h(labindex);      % unique step size corresponding to worker index. 
    t=t0:ht:tf;
    y=zeros(length(t),1);
    z=zeros(length(t),1);
    y(1)=y0;
    z(1)=y0;
    
    for i=1:length(t);
        % Integrate one step.
        z(i+1)=z(i)+ht*(lambda*z(i)+z(i)^p);     % FE.
    
        y1=y(i)+ht*lambda*y(i);      % SpFE.
        y(i+1)=blowup(y1,ht,p);
    end
    % calculate absolute errors.
    ge=abs(y(i)-yexact);
    geFE=abs(z(i)-yexact);
    GE=[ge,geFE];



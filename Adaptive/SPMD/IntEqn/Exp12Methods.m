function [GE, GESp]=Exp12Methods(atol,rtol,t0,tf,h0,y0,n,hx,delta,aij,yexact)
format long;
ge=0;
geSp=0;
options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,...
            'Display','off');
if (labindex<=4)
    k=1;
    t(1)=t0;
    ht=h0;
    hfailed=false;
    frac=0.9;
    FE(:,1)=y0;
    EMR(:,1)=y0;
    while (abs(t(k)-tf) >= 1e-16)
        t(k);
        FE(:,k+1)=FE(:,k)+ht*f(k*ht,FE(:,k),n,...
                    hx,delta,aij);
    
        % explicit midpoint rule.
        ztmp=EMR(:,k)+0.5*ht*f(k*ht,EMR(:,k),n,hx,delta,aij);
        ztmp(1)=0;
        ztmp(n)=0;
        
        EMR(:,k+1)=EMR(:,k)+ht*f(k*ht,ztmp,n,hx,delta,aij);
        EMR(1,k+1)=0;
        EMR(n,k+1)=0;
        
        etol=atol(labindex)+rtol(labindex)*norm(FE(:,k+1),inf);  % tolerance.
        epsilon=norm(FE(:,k+1)-EMR(:,k+1),inf);  % local error estimate.
    
        if (epsilon<=frac*etol) % accepted local step.
            FE(:,k+1)=EMR(:,k+1);   % set lower order sol'n to higher order sol'n.
            t(k+1)=t(k)+ht;
            s =(etol/epsilon)^(0.5);
            ht= s*ht;   % optimal step size
            k=k+1;
            if t(k)+ht>tf
                ht=tf-t(k);
            end
        else                  % rejected local step.
            hfailed=true;
            s = (etol/epsilon)^(0.5);  %scalar for new h
            htmp= s*ht; % optimal step size
            if hfailed % do  not allow step to increase if last step failed
                if (htmp>=ht)
                    ht=0.9*ht;
                else
                    ht=htmp;
                end
            end
        end
    end
    ge=norm(yexact-FE(:,k),inf);
elseif (labindex>4 && labindex<=8)
    k=1;
    t(1)=t0;
    ht=h0;
    hfailed=false;
    frac=0.9;
    SpFE(:,1)=y0;
    Str(:,1)=y0;
    while (abs(t(k)-tf) >= 1e-16)
        % one integration step. 
        % B-method.
        ztmp=zeros(n,1);       % SpFE.
        ztmp(2:n-1)=-log(exp(-SpFE(2:n-1,k))-delta*ht);
        SpFE(2:n-1,k+1)=ztmp(2:n-1)+ht*(ztmp(1:n-2)-2*...
            ztmp(2:n-1)+ztmp(3:n))/(hx^2);
  
        % Strang.
        y15=zeros(n,1);
        y15(2:n-1)=-log(exp(-Str(2:n-1,k))-delta*ht/2);
        y25=ExpSecOrd(y15,ht/2,n,hx);
        y35=intdiff(y25,ht,hx,n,aij);
        y45=ExpSecOrd(y35,ht/2,n,hx);
        Str(2:n-1,k+1)=-log(exp(-y45(2:n-1))-delta*ht/2);
    
        etol=atol(labindex-4)+rtol(labindex-4)*norm(SpFE(:,k+1),inf);  % tolerance.
        epsilon=norm(SpFE(:,k+1)-Str(:,k+1),inf);  % local error estimate.
    
        if (epsilon<=frac*etol) % accepted local step.
            SpFE(:,k+1)=Str(:,k+1);   % set lower order sol'n to higher order sol'n.
            t(k+1)=t(k)+ht;
            s =(etol/epsilon)^(0.5);
            ht= s*ht;   % optimal step size
            k=k+1;
            if t(k)+ht>tf
                ht=tf-t(k);
            end
        else                  % rejected local step.
            hfailed=true;
            s = (etol/epsilon)^(0.5);  %scalar for new h
            htmp= s*ht; % optimal step size
            if hfailed % do  not allow step to increase if last step failed
                if (htmp>=ht)
                    ht=0.9*ht;
                else
                    ht=htmp;
                end
            end
        end
    end
    geSp=norm(yexact-SpFE(:,k),inf);
end
codist=codistributor1d(1,[1 1 1 1 1 1 1 1],[8 1]);
GE=codistributed.build(ge,codist);
GESp=codistributed.build(geSp,codist);
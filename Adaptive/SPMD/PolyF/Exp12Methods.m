function [GE, GESp]=Exp12Methods(atol,rtol,t0,tf,h0,y0,n,hx,delta,alpha,p,yexact)
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
        % one integration step.  
        FE(:,k+1)=FE(:,k)+ht*fbeck(k*ht,...         % FE.
            FE(:,k),n,hx,delta,alpha,p);
    
        ztmp=EMR(:,k)+0.5*ht*fbeck(k*ht,...   %EMR.
            EMR(:,k),n,hx,delta,alpha,p);
        ztmp(1)=0;
        ztmp(n)=0;
        EMR(:,k+1)=EMR(:,k)+ht*fbeck(k*ht,ztmp,n,hx,delta,alpha,p);
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
    SoSpFE(:,1)=y0;
    while (abs(t(k)-tf) >= 1e-16)
        % one integration step. 
        ztmp=zeros(n,1);       % SpFE.
        ztmp(2:n-1)=((SpFE(2:n-1,k)+alpha*ones(1,n-2)').^(1-p)+(1-p)*delta*...
            ht).^(1/(1-p))-alpha;
        SpFE(2:n-1,k+1)=ztmp(2:n-1)+ht*(ztmp(1:n-2)-2*ztmp(2:n-1)+ztmp(3:n))...
            /(hx^2);
  
        v=fsolve(@(x) SoPsi(x,SoSpFE(:,k),ht,n,hx,...       % SpSoFE.
            delta,alpha,p),SoSpFE(:,k),options);
        SoSpFE(2:n-1,k+1)=((v(2:n-1)+(0.5*ht/hx^2)*(v(1:n-2)-...
            2*v(2:n-1)+v(3:n))+alpha).^(1-p)+0.5*(1-p)*delta*ht)...
            .^(1/(1-p))-alpha;
    
        etol=atol(labindex-4)+rtol(labindex-4)*norm(SpFE(:,k+1),inf);  % tolerance.
        epsilon=norm(SpFE(:,k+1)-SoSpFE(:,k+1),inf);  % local error estimate.
    
        if (epsilon<=frac*etol) % accepted local step.
            SpFE(:,k+1)=SoSpFE(:,k+1);   % set lower order sol'n to higher order sol'n.
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

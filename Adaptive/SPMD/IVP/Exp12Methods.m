function [GE, GESp]=Exp12Methods(atol,rtol,t0,tf,h0,y0,lambda,p)
ge=0;
geSp=0;
lamp=lambda*(1-p);
if (labindex<=4)
    i=1;
    t(1)=t0;
    ht=h0;
    hfailed=false;
    frac=0.9;
    FE(1)=y0;
    EMR(1)=y0;
    while (abs(t(i)-tf) >= 1e-16)
        % one integration step.  
        FE(i+1)=FE(i)+ht*(lambda*FE(i)+FE(i)^p);
    
        ztmp=EMR(i)+0.5*ht*(lambda*EMR(i)+EMR(i)^p);             
        EMR(i+1)=EMR(i)+ht*(lambda*ztmp+ztmp^p);
    
        etol=atol(labindex)+rtol(labindex)*norm(FE(i+1));  % tolerance.
        epsilon=FE(i+1)-EMR(i+1);  % local error estimate.
    
        if (abs(epsilon)<=frac*etol) % accepted local step.
            FE(i+1)=EMR(i+1);   % set lower order sol'n to higher order sol'n.
            t(i+1)=t(i)+ht;
            s =(etol/(norm(epsilon)))^(1/2);
            ht= s*ht;   % optimal step size
            i=i+1;
            if t(i)+ht>tf
                ht=tf-t(i);
            end
        else                  % rejected local step.
            hfailed=true;
            s = (etol/(norm(epsilon)))^(1/2);  %scalar for new h
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
    yexact=(((y0^(1-p)+(1/lambda))*exp(lamp*t(i)))-(1/lambda))^(1/(1-p));
    ge=abs(yexact-FE(i));
elseif (labindex>4 && labindex<=8)
    i=1;
    t(1)=t0;
    ht=h0;
    hfailed=false;
    frac=0.9;
    SpFE(1)=y0;
    SoSpFE(1)=y0;
    while (abs(t(i)-tf) >= 1e-16)
        % one integration step. 
        ytmp=SpFE(i)+ht*lambda*SpFE(i);   
        SpFE(i+1)=((1-p)*ht+ytmp^(1-p))^(1/(1-p));
    
        SoSpFE(i+1)=((1-p)*0.5*ht+(((1+0.5*lambda*ht)/(1-0.5*lambda*ht))...
            ^(1-p))*(0.5*(1-p)*ht+SoSpFE(i)^(1-p)))^(1/(1-p));
    
        etol=atol(labindex-4)+rtol(labindex-4)*norm(SpFE(i+1));  % tolerance.
        epsilon=SpFE(i+1)-SoSpFE(i+1);  % local error estimate.
    
        if (abs(epsilon)<=frac*etol) % accepted local step.
            SpFE(i+1)=SoSpFE(i+1);   % set lower order sol'n to higher order sol'n.
            t(i+1)=t(i)+ht;
            s =(etol/(norm(epsilon)))^(1/2);
            ht= s*ht;   % optimal step size
            i=i+1;
            if t(i)+ht>tf
                ht=tf-t(i);
            end
        else                  % rejected local step.
            hfailed=true;
            s = (etol/(norm(epsilon)))^(1/2);  %scalar for new h
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
    yexact=(((y0^(1-p)+(1/lambda))*exp(lamp*t(i)))-(1/lambda))^(1/(1-p));
    geSp=abs(yexact-SpFE(i));
end
codist=codistributor1d(1,[1 1 1 1 1 1 1 1],[8 1]);
GE=codistributed.build(ge,codist);
GESp=codistributed.build(geSp,codist);

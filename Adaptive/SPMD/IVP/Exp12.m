clear all;
format long;
tic

p=4;
lambda=-1;
lamp=lambda*(1-p);

y0=2;            % initial condition.
FE(1)=y0;
EMR(1)=y0;

Tb=(1/lamp)*log(1/(lambda*y0^(1-p)+1));  % theoretical blow-up time.
tf=0.999*Tb;
rtol=10^-7;
atol=10^-7;
ht=0.0006;
frac=0.9;
i=1;
hfailed=false;
t(1)=0;

while (abs(t(i)-tf) >= 1e-16)
    % one integration step.  
    FE(i+1)=FE(i)+ht*(lambda*FE(i)+FE(i)^p);
    
    ztmp=EMR(i)+0.5*ht*(lambda*EMR(i)+EMR(i)^p);             
    EMR(i+1)=EMR(i)+ht*(lambda*ztmp+ztmp^p);
    
    etol=atol+rtol*norm(FE(i+1));  % tolerance.
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
        ge=abs(yexact-FE(i))
        t(i)/Tb


toc
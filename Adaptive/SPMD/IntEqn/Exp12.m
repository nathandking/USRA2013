clear all;
format long;
tic

t0=0;                    % initial time.
Tb=0.1901947;
tf=0.999*Tb;
n=25;  % number of spatial discretization nodes.
delta=3;
ht=0.0001;  % initial time step size.

opt=odeset('AbsTol',1e-12,'RelTol',1e-12); 
x=linspace(-1,1,n);    % define spatial nodes on the inteval [-1,1].
aij=a(x,x);
hx=2/(n-1);            % uniform distance between spatial nodes.
y0=cos(pi*x/2);           % initial condition.
sol=ode45(@(t,y) f(t,y,...   % compute reference solution for each n.
    n,hx,delta,aij),[t0 tf],y0,opt); 
yexact=deval(sol,tf);
    
FE(:,1)=y0;
EMR(:,1)=y0;
rtol=10^-5;
atol=10^-5;
frac=0.9;
k=1;
hfailed=false;
t(1)=t0;

while (abs(t(k)-tf) >= 1e-16)
    % one integration step.  
    FE(:,k+1)=FE(:,k)+ht*f(k*ht,FE(:,k),n,...
                    hx,delta,aij);
    
    % explicit midpoint rule.
    ztmp=EMR(:,k)+0.5*ht*f(k*ht,EMR(:,k),n,hx,delta,aij);
    ztmp(1)=0;
    ztmp(n)=0;
        
    EMR(:,k+1)=EMR(:,k)+ht*f(k*ht,ztmp,n,hx,delta,aij);
    EMR(1,k+1)=0;
    EMR(n,k+1)=0;
    
    etol=atol+rtol*norm(FE(:,k+1),inf);  % tolerance.
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
ge=norm(yexact-FE(:,k),inf)
fract=t(k)/Tb;
toc
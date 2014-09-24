clear all;
format long;
tic

t0=0;                    % initial time.
Tb=0.1661945;
tf=0.99*Tb;
n=25;  % number of spatial discretization nodes.
delta=3;
ht=0.0001;  % initial time step size.

opt=odeset('AbsTol',1e-12,'RelTol',1e-12); 
x=linspace(-1,1,n);    % define spatial nodes on the inteval [-1,1].          
hx=2/(n-1);            % uniform distance between spatial nodes.
y0=cos(pi*x/2);           % initial condition.
sol=ode45(@(t,y) fbeck(t,y,...   % compute reference solution for each n.
    n,hx,delta),[t0 tf],y0,opt); 
yexact=deval(sol,tf);
    
SpFE(:,1)=y0;
SoSpFE(:,1)=y0;
rtol=10^-5;
atol=10^-5;
frac=0.9;
k=1;
hfailed=false;
t(1)=t0;
options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,...
            'Display','off');

while (abs(t(k)-tf) >= 1e-16)
    % one integration step.  
    ztmp=zeros(n,1);       % SpFE.
    ztmp(2:n-1)=-log(exp(-SpFE(2:n-1,k))-delta*ht);
    SpFE(2:n-1,k+1)=ztmp(2:n-1)+ht*(ztmp(1:n-2)-2*...
                    ztmp(2:n-1)+ztmp(3:n))/(hx^2);
  
    v=fsolve(@(x) SoPsi(x,SoSpFE(:,k),ht,n,delta,hx),SoSpFE(:,k),...
        options);
    SoSpFE(2:n-1,k+1)=-log(exp(-v(2:n-1)-(ht/(2*hx^2))*...
        (v(1:n-2)-2*v(2:n-1)+v(3:n)))-0.5*delta*ht);
    
    etol=atol+rtol*norm(SpFE(:,k+1),2);  % tolerance.
    epsilon=norm(SpFE(:,k+1)-SoSpFE(:,k+1),2);  % local error estimate.
    
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
ge=norm(yexact-SpFE(:,k),2)
t(k)/Tb


toc
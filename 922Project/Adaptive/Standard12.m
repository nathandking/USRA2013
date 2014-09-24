% This is the standard adaptive time stepping method BE-CN.

clear all;
format long;
tic

t(1)=0;                    % initial time.
Tb=0.494;
tf=0.99*Tb;
n=51;  % number of spatial discretization nodes.
alpha=1;
k=0.0001;  % initial time step size.

options=optimset('TolFun',1e-10);
opt=odeset('AbsTol',1e-12,'RelTol',1e-12); 
x=linspace(-1,1,n);    % define spatial nodes on the interval [-1,1].          
hx=2/(n-1);            % uniform distance between spatial nodes.
beta=4;
y0=beta*cos(pi*x/2);           % initial condition.
sol=ode15s(@(t,y) f(t,y,hx,alpha),[t(1) tf],y0,opt);   % compute reference solution for each n.
    
yBE(:,1)=y0;
yCN(:,1)=y0;
rtol=10^-6;
atol=10^-6;
frac=0.9;
i=1;
hfailed=false;
ge=0;

while (ge <= 1e-5)
    % one integration step.  
    
    % baciward Euler.
    yBE(:,i+1)=fsolve(@(x) BE(x,yBE(:,i),k,hx,alpha),yBE(:,i),options);
    
    % Crank-Nicolson.
    yCN(:,i+1)=fsolve(@(x) CN(x,yCN(:,i),k,hx,alpha),yCN(:,i),options);
    
    etol=atol+rtol*norm(yBE(:,i+1),inf);  % tolerance.
    epsilon=norm(yBE(:,i+1)-yCN(:,i+1),inf);  % local error estimate.
    
    if (epsilon<=frac*etol) % accepted local step.
        yBE(:,i+1)=yCN(:,i+1);   % set lower order sol'n to higher order sol'n.
        t(i+1)=t(i)+k;
        s =(etol/epsilon)^(0.5);
        k= s*k;   % optimal step size
        yexact=deval(sol,t(i+1));
        ge(i)=norm(yexact-yBE(:,i+1),inf);
        i=i+1;
    else                  % rejected local step.
        hfailed=true;
        s = (etol/epsilon)^(0.5);  %scalar for new h
        kmp= s*k; % optimal step size
        if hfailed % do  not allow step to increase if last step failed
            if (kmp>=k)
                k=0.9*k;
            else
                k=kmp;
            end
        end
    end
end        
t(i)/Tb


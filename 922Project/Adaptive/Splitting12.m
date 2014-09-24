% This is the splitting adaptive time stepping method SpBE-SpCN.

clear all;
format long;
tic

tSp(1)=0;                    % initial time.
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
sol=ode15s(@(t,y) f(t,y,hx,alpha),[tSp(1) tf],y0,opt);   % compute reference solution for each n.
    
ySpBE(:,1)=y0;
vSpBE(:,1)=y0;
rtol=10^-6;
atol=10^-6;
frac=0.9;
i=1;
hfailed=false;
geSp=0;

while (geSp <= 1e-5)
    % one integration step.  
    
    % First-Order Splitting.
    y1=fsolve(@(x) SpBE(x,ySpBE(:,i),k,hx),ySpBE(:,i),options);
        
    ySpBE(1,i+1)=0;
    ySpBE(2:n-1,i+1)=(y1(2:n-1).^(-alpha)-alpha*k).^(1/-alpha);
    ySpBE(n,i+1)=0;
    
    % Second-Order Splitting.
    v1=zeros(n,1);
    v1(2:n-1)=(vSpBE(2:n-1,i).^(-alpha)-0.5*alpha*k).^(1/-alpha);
        
    v2=zeros(n,1);
    v2(2:n-1)=v1(2:n-1)+0.5*k*(v1(1:n-2)-2*v1(2:n-1)+v1(3:n))/(hx^2);
       
    v3=fsolve(@(x) SpBE(x,v2,0.5*k,hx),v2,options);
        
    vSpBE(1,i+1)=0;
    vSpBE(2:n-1,i+1)=(v3(2:n-1).^(-alpha)-0.5*alpha*k).^(1/-alpha);
    vSpBE(n,i+1)=0;
    
    etol=atol+rtol*norm(ySpBE(:,i+1),inf);  % tolerance.
    epsilon=norm(ySpBE(:,i+1)-vSpBE(:,i+1),inf);  % local error estimate.
    
    if (epsilon<=frac*etol) % accepted local step.
        ySpBE(:,i+1)=vSpBE(:,i+1);   % set lower order sol'n to higher order sol'n.
        tSp(i+1)=tSp(i)+k;
        s =(etol/epsilon)^(0.5);
        k= s*k;   % optimal step size
        yexact=deval(sol,tSp(i+1));
        geSp=norm(yexact-ySpBE(:,i+1),inf)
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
tSp(i)/Tb
figure(1)
plot(1:i-1,geSp,'ro')
xlabel('$n$','Interpreter','latex','FontSize',16)
ylabel('Global Error','Interpreter','latex','FontSize',16)
figure(2)
plot(1:i-1,diff(tSp),'o')
xlabel('$n$','Interpreter','latex','FontSize',16)
ylabel('Step Size','Interpreter','latex','FontSize',16)


toc
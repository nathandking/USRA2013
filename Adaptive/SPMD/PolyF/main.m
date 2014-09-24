clear all;
format long;

if matlabpool('size')==0          % if no matlab pool exists.
    matlabpool('open',8);
elseif matlabpool('size')~=8     % if matlab pool does not have 11 cores.
    matlabpool('close');
    matlabpool('open',8);
end

n=51;                                    % spatial discretization nodes.
p=2;
delta=3;
alpha=2;
t0=0;                                    % initial time.
h0=0.0001;                               % initial step size.
Tb=0.1209365;                                     % numerical blow-up time.
tf=[0.5*Tb,0.95*Tb,0.9999*Tb];           % final time of integration.
rtol=[1e-1; 1e-3; 1e-5; 1e-7];           % relative error tolerances.
atol=[1e-1; 1e-3; 1e-5; 1e-7];           % absolute error tolerances.
ge=zeros(3,4);
geSp=zeros(3,4);

opt=odeset('AbsTol',1e-12,'RelTol',1e-12); 
x=linspace(-1,1,n);    % define spatial nodes on the interval [-1,1].          
hx=2/(n-1);            % uniform distance between spatial nodes.
y0=cos(pi*x/2);           % initial condition.
sol=ode45(@(t,y) fbeck(t,y,...   % compute reference solution for each n.
    n,hx,delta,alpha,p),[t0 tf],y0,opt); 
yexact=deval(sol,tf);

for i=1:3
spmd
  [GE,GESp]=Exp12Methods(atol,rtol,t0,tf(i),h0,y0,n,hx,delta,alpha,p,yexact(:,i)); 
end

GE=gather(GE);
GESp=gather(GESp);
ge(i,1:4)=GE(1:4);
geSp(i,1:4)=GESp(5:8);
end

for i=1:4
fprintf('%1.0e & %10.4e & %10.4e & %10.4e & %10.4e & %10.4e & %10.4e\\\\\n'...
    ,rtol(i),ge(1,i),geSp(1,i),ge(2,i),geSp(2,i),ge(3,i),geSp(3,i));
end
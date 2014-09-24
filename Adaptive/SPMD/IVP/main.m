clear all
format long

if matlabpool('size')==0          % if no matlab pool exists.
    matlabpool('open',8);
elseif matlabpool('size')~=8     % if matlab pool does not have 11 cores.
    matlabpool('close');
    matlabpool('open',8);
end

p=4;
lambda=-1;
lamp=lambda*(1-p);
y0=2;                                    % initial condition.
t0=0;                                    % initial time.
h0=0.0001;                               % initial step size.
Tb=(1/lamp)*log(1/(lambda*y0^(1-p)+1));  % theoretical blow-up time.
tf=[0.5*Tb,0.95*Tb,0.9999*Tb];           % final time of integration.
rtol=[1e-1; 1e-3; 1e-5; 1e-7];           % relative error tolerances.
atol=[1e-1; 1e-3; 1e-5; 1e-7];           % absolute error tolerances.
ge=zeros(3,4);
geSp=zeros(3,4);

for i=1:3
spmd
  [GE,GESp]=Exp12Methods(atol,rtol,t0,tf(i),h0,y0,lambda,p); 
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
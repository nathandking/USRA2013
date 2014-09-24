clear all;
format long;
%-------------------------------------------------------------------------
%   Implementation of a second order B-method detailed by Beck in her
%   thesis. This B-method, of splitting type, uses forward Euler for 
%   Phi^[2]. The PDE solved is given by
%             u_t=u_xx+delta*F(u),
%   with Dirichlet boundary conditions and F(u)=exp(u).
%-------------------------------------------------------------------------

global n delta hx;

t(1)=0;                % initial time.
tf=0.1;             % final time of integration.
n=301;                  % number of spatial nodes.
x=linspace(-1,1,n);    % define spatial nodes on the inteval [-1,1].          
hx=2/(n-1);            % uniform distance between spatial nodes.
y(:,1)=cos(pi*x/2);    % initial condition.
delta=3;
ht=[0.0002 0.000125 0.0001 0.00005 0.000025];  % time step sizes.

options=optimset('TolFun',1e-16);
opt=odeset('AbsTol',1e-12,'RelTol',1e-12);
[T yexact]=ode45(@f,[t(1) tf],y(:,1),opt);  % reference solution.
  
% complete integration using different time step sizes.
for j=1:5
    t=0:ht(j):tf;
    for i=1:length(t)
        % Integrate one step.
        y15(1,i+1)=0;
        y15(2:n-1,i+1)=-log(exp(-y(2:n-1,i))-delta*ht(j)/2);
        y15(n,i+1)=0;
        
%         y25(1,i+1)=0;
%         y25(2:n-1,i+1)=y15(2:n-1,i+1)+(ht(j)/2)*(y15(1:n-2,i+1)-...
%             2*y15(2:n-1,i+1)+y15(3:n,i+1))/(hx^2);
%         y25(n,i+1)=0;

        y25(:,i+1)=ImpSecOrd(y15(:,i+1),ht(j)/2);
        
        y35(:,i+1)=intdiff(y25(:,i+1),ht(j));
        
%         y45(:,i+1)=fsolve(@(x) SoPsi(x,y35(:,i+1),ht(j)/2),y35(:,i+1),options);
        y45(:,i+1)=ImpSecOrd(y35(:,i+1),ht(j)/2);
       
        y(1,i+1)=0;
        y(2:n-1,i+1)=-log(exp(-y45(2:n-1,i+1))-delta*ht(j)/2);
        y(n,i+1)=0;
    end
    ge(j)=norm(y(:,i)-yexact(end,:)',inf);
end

% print results.
ge'
fprintf('Order test:\n')
fprintf('%10.8f\n',ge(3)/ge(4))

figure(1)
plot(x,y(:,length(t)),'.',x,yexact(end,:))
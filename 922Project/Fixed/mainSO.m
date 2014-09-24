clear all;
format long;
%-------------------------------------------------------------------------
%   Implementation of a second order splitting method detailed by Beck in 
%   her thesis. This splitting method uses backward Euler for 
%   Phi^*[2]. The PDE solved is given by
%             u_t=u_xx+F(u),
%   with Dirichlet boundary conditions and F(u) = u^{1+alpha}.
%-------------------------------------------------------------------------

global n hx alpha;

t(1)=0;                % initial time.
tf=0.493;               % final time of integration. (Tb=4.944521e-01)
n=51;                  % number of spatial nodes.
x=linspace(-1,1,n);    % define spatial nodes on the inteval [-1,1].          
hx=2/(n-1);            % uniform distance between spatial nodes.
alpha=1;
beta=4;

ySpBE(:,1)=beta*cos(0.5*pi*x);    % initial condition.
yCN(:,1)=ySpBE(:,1);
% yStrang(:,1)=ySpBE(:,1);

k=[0.0001 0.00005 0.000025 0.0000125];  % time step sizes.

options=optimset('TolFun',1e-10);
opt=odeset('AbsTol',1e-12,'RelTol',1e-12);
[T yexact]=ode15s(@Ff,[t(1) tf],ySpBE(:,1),opt);  % reference solution.

% complete integration using different time step sizes.
for j=1:4
    t=0:k(j):tf;
    for i=1:length(t)
        % integrate one step.
        
        % Crank-Nicolson.
        yCN(:,i+1)=fsolve(@(x) FCN(x,yCN(:,i),k(j)),yCN(:,i),options);
        
        % Second-Order Splitting.
        y1=zeros(n,1);
        y1(2:n-1)=(ySpBE(2:n-1,i).^(-alpha)-0.5*alpha*k(j)).^(1/-alpha);
        
        y2=zeros(n,1);
        y2(2:n-1)=y1(2:n-1)+0.5*k(j)*(y1(1:n-2)-2*y1(2:n-1)+y1(3:n))/(hx^2);
       
        y3=fsolve(@(x) FSpBE(x,y2,0.5*k(j)),y2,options);
        
        ySpBE(1,i+1)=0;
        ySpBE(2:n-1,i+1)=(y3(2:n-1).^(-alpha)-0.5*alpha*k(j)).^(1/-alpha);
        ySpBE(n,i+1)=0;
        
        % Strang Splitting (same as other B-method because using CN).
%         v1=zeros(n,1);
%         v1(2:n-1)=(ySpBE(2:n-1,i).^(-alpha)-0.5*alpha*k(j)).^(1/-alpha);
%         
%         v2=fsolve(@(x) CNStrang(x,v1,k(j)),v1,options);
%         
%         yStrang(1,i+1)=0;
%         yStrang(2:n-1,i+1)=(v2(2:n-1).^(-alpha)-0.5*alpha*k(j)).^(1/-alpha);
%         yStrang(n,i+1)=0;
        
    end
    ge(j)=norm(ySpBE(:,i)-yexact(end,:)',inf);
    geCN(j)=norm(yCN(:,i)-yexact(end,:)',inf);
%     geStrang(j)=norm(yStrang(:,i)-yexact(end,:)',inf);

    if j==1
        clf;
        figure(1)
        plot(t(4301:end),yCN(26,4301:end-1),'bd')
        hold on
        plot(t(4301:end),ySpBE(26,4301:end-1),'ro')
        hold on
        plot(T(1235:end),yexact(1235:end,26),'k')
        xlabel('$t$','Interpreter','latex','FontSize',16)
        title('Solution at $x=0$','Interpreter','latex','FontSize',16)
        axis([0.485 0.4935 100 800])
        l=legend('CN','SpCN','  ode15s',2);
        set(l,'Interpreter','latex','FontSize',16);
        set(l,'position',[0.15 0.78 0.2 0.1])
    end
end

% print results.
f1=figure(2);         % create convergence plots.
loglog(k,geCN,'-bd',k,ge,'-ro');
xlabel('$k$','Interpreter','latex','FontSize',16);
ylabel('Error','Interpreter','latex','FontSize',16);
axis([1e-5 1.2e-4 1e-5 1e1]);
print(f1,'-deps','~/plot.eps');
l=legend('CN','SpCN',2);
set(l,'Interpreter','latex','FontSize',14);
set(l,'position',[0.15 0.8 0.2 0.1])
hold on
plot([2e-5 9e-5],[0.00256 0.05184],'k--')
h=text(1.35e-5,1.2e-3,'Slope=2','FontSize',14,'Interpreter','latex');
set(h,'rotation',17)



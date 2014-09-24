% load ExpFFO.txt 
% h=[0.00005, 0.000025, 0.0000125, 0.000008, 0.000005];
% figure(1)
% loglog(h,ExpFFO(1,:),'-sk',h,ExpFFO(2,:),'-*k',h,ExpFFO(3,:),'-ok',h,ExpFFO(4,:),'-dk')
% xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
% ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
% legend('FE','FEA','SpFEA','SpFE',2);
% axis([4e-6 6.5e-5 4e-4 1e0]);
% hold on
% plot([1e-5 4.5e-5],[0.01 0.045],'k--')
% h=text(6e-6,6.2e-3,'Slope=1','FontSize',16,'Interpreter','latex');
% set(h,'rotation',15)
 
 
load ExpFSO.txt
h=[0.0002; 0.000125; 0.0001; 0.00005; 0.000025];
figure(2)
loglog(h,ExpFSO(1,:),'-sk',h,ExpFSO(2,:),'-*k',h,ExpFSO(3,:),'-ok',h,ExpFSO(4,:),'-dk')
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2);
axis([2e-5 2.5e-4 5e-6 1e-1]);
hold on
plot([5e-5 1.8e-4],[5e-4 0.0065],'k--')
h=text(3.2e-5,2.1e-4,'Slope=2','FontSize',16,'Interpreter','latex');
set(h,'rotation',22)

t0=0;                    % initial time.
tf=0.1660;               % final time of integration.
n=51;
delta=3;
opt=odeset('AbsTol',1e-12,'RelTol',1e-12);
x=linspace(-1,1,n);    % define spatial nodes on the inteval [-1,1].          
hx=2/(n-1);            % uniform distance between spatial nodes.
y0=cos(pi*x/2);           % initial condition.
sol=ode45(@(t,y) fbeck(t,y,n,hx,delta),[t0 tf],y0,opt); 
yexact=deval(sol,tf);

% load ExpFFO.txt 
% ExpFFO=ExpFFO./norm(yexact,inf);
% h=[0.00005, 0.000025, 0.0000125, 0.000008, 0.000005];
% figure(3)
% loglog(h,ExpFFO(1,:),'-sk',h,ExpFFO(2,:),'-*k',h,ExpFFO(3,:),'-ok',h,ExpFFO(4,:),'-dk')
% xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
% ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
% legend('FE','FEA','SpFEA','SpFE',2);
% axis([4e-6 6.5e-5 5e-5 2e-1]);
% hold on
% plot([1e-5 4.5e-5],[0.001 0.0045],'k--')
% h=text(6e-6,6.2e-4,'Slope=1','FontSize',16,'Interpreter','latex');
% set(h,'rotation',15)
 
 
load ExpFSO.txt
ExpFSO=ExpFSO./norm(yexact,inf);
h=[0.0002; 0.000125; 0.0001; 0.00005; 0.000025];
figure(4)
loglog(h,ExpFSO(1,:),'-sk',h,ExpFSO(2,:),'-*k',h,ExpFSO(3,:),'-ok',h,ExpFSO(4,:),'-dk')
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2);
axis([2e-5 2.5e-4 5e-7 1.5e-2]);
hold on
plot([5e-5 1.8e-4],[5e-5 0.00065],'k--')
h=text(3.2e-5,2.1e-5,'Slope=2','FontSize',16,'Interpreter','latex');
set(h,'rotation',22)
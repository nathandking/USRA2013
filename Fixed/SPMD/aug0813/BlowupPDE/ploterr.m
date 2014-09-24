% load BlowupPDEFOSp51.txt 
% h=[0.0001; 0.00005; 0.000025; 0.0000125; 0.000008];
% figure(1)
% loglog(h,BlowupPDEFOSp51(1,:),'-sk',h,BlowupPDEFOSp51(2,:),'-*k',h,BlowupPDEFOSp51(3,:),'-ok',h,BlowupPDEFOSp51(4,:),'-dk')
% xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
% ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
% legend('FE','FEA','SpFEA','SpFE',2);
% axis([5e-6 2e-4 1e-3 1e1]);
% hold on
% plot([2e-5 8e-5],[0.08 0.32],'k--')
% h=text(1e-5,4e-2,'Slope=1','FontSize',16,'Interpreter','latex');
% set(h,'rotation',17)
 
 
load BlowupPDESOSp51.txt
h=[0.0002; 0.0001; 0.00005; 0.000025; 0.0000125];
figure(2)
loglog(h,BlowupPDESOSp51(1,:),'-sk',h,BlowupPDESOSp51(2,:),'-*k',h,BlowupPDESOSp51(3,:),'-ok',h,BlowupPDESOSp51(4,:),'-dk')
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e-1]);
hold on
plot([3e-5 1.5e-4],[4.5e-5 1.1e-3],'k--')
h=text(1.7e-5,1.6e-5,'Slope=2','FontSize',16,'Interpreter','latex');
set(h,'rotation',19)

t0=0;                    % initial time.
tf=0.1150;               % final time of integration.
n=51;  % number of spatial discretization nodes.
delta=3;
p=2;
alpha=2;
opt=odeset('AbsTol',1e-12,'RelTol',1e-12);
    x=linspace(-1,1,n);    % define spatial nodes on the inteval [-1,1].          
    hx=2/(n-1);            % uniform distance between spatial nodes.
    y0=cos(pi*x/2);           % initial condition.
    sol=ode45(@(t,y) fbeck(t,y,n,hx,delta,alpha,p),[t0 tf],y0,opt); 
    yexact=deval(sol,tf);

% load BlowupPDEFOSp51.txt 
% BlowupPDEFOSp51=BlowupPDEFOSp51./norm(yexact,inf);
% h=[0.0001; 0.00005; 0.000025; 0.0000125; 0.000008];
% figure(3)
% loglog(h,BlowupPDEFOSp51(1,:),'-sk',h,BlowupPDEFOSp51(2,:),'-*k',h,BlowupPDEFOSp51(3,:),'-ok',h,BlowupPDEFOSp51(4,:),'-dk')
% xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
% ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
% legend('FE','FEA','SpFEA','SpFE',2);
% axis([5e-6 2e-4 1e-5 1e-1]);
% hold on
% plot([2e-5 8e-5],[0.0008 0.0032],'k--')
% h=text(1e-5,4e-4,'Slope=1','FontSize',16,'Interpreter','latex');
% set(h,'rotation',17)
 
 
load BlowupPDESOSp51.txt
BlowupPDESOSp51=BlowupPDESOSp51./norm(yexact,inf);
h=[0.0002; 0.0001; 0.00005; 0.000025; 0.0000125];
figure(4)
loglog(h,BlowupPDESOSp51(1,:),'-sk',h,BlowupPDESOSp51(2,:),'-*k',h,BlowupPDESOSp51(3,:),'-ok',h,BlowupPDESOSp51(4,:),'-dk')
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2);
axis([1e-5 2.5e-4 3e-9 1.5e-3]);
hold on
plot([3e-5 1.5e-4],[4.5*10^(-6.8) 1.1*10^(-4.8)],'k--')
h=text(1.7e-5,1.6*10^(-6.8),'Slope=2','FontSize',16,'Interpreter','latex');
set(h,'rotation',19)
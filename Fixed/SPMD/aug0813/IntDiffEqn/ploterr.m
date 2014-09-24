load IntFO.txt 
h=[0.0001, 0.00005, 0.000025, 0.0000125, 0.000008];
figure(1)
loglog(h,IntFO(1,:),'-sk',h,IntFO(2,:),'-*k',h,IntFO(3,:),'-ok',h,IntFO(4,:),'-dk')
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('FE','FEA','SpFEA','SpFE',2);
axis([6e-6 1.5e-4 1e-3 2e0]);
hold on
plot([2e-5 8e-5],[0.04 0.16],'k--')
h=text(1.1e-5,2.2e-2,'Slope=1','FontSize',16,'Interpreter','latex');
set(h,'rotation',19)
 
 
% load IntSO.txt
% h=[0.0002; 0.000125; 0.0001; 0.00005; 0.000025];
% figure(2)
% loglog(h,IntSO(1,:),'-sk',h,IntSO(2,:),'-*k',h,IntSO(3,:),'-ok')
% xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
% ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
% legend('EMR','TR','Strang',2);
% axis([2e-5 2.5e-4 2e-6 8e-2]);
% hold on
% plot([4.5e-5 1.8e-4],[0.00015188 0.0024],'k--')
% h=text(2.9e-5,6.7e-5,'Slope=2','FontSize',16,'Interpreter','latex');
% set(h,'rotation',20)

t0=0;                  % initial time.
tf=0.19;             % final time of integration.
n=51;                  % number of spatial nodes.
delta=3;
% compute the reference solution.
opt=odeset('AbsTol',1e-12,'RelTol',1e-12);
    x=linspace(-1,1,n);    % define spatial nodes on the inteval [-1,1].          
    hx=2/(n-1);            % uniform distance between spatial nodes.
    y0=cos(pi*x/2);           % initial condition.
    aij=a(x,x);
    sol=ode45(@(t,y) f(t,y,n,hx,delta,aij),[t0 tf],y0,opt); 
    yexact=deval(sol,tf);

load IntFO.txt 
IntFO=IntFO./norm(yexact,inf);
h=[0.0001, 0.00005, 0.000025, 0.0000125, 0.000008];
figure(3)
loglog(h,IntFO(1,:),'-sk',h,IntFO(2,:),'-*k',h,IntFO(3,:),'-ok',h,IntFO(4,:),'-dk')
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('FE','FEA','SpFEA','SpFE',2);
axis([6e-6 1.5e-4 1*10^(-3.8) 2*10^(-0.8)]);
hold on
plot([2e-5 8e-5],[0.004*10^(0.2) 0.016*10^(0.2)],'k--')
h=text(1.1e-5,2.2*10^(-2.8),'Slope=1','FontSize',16,'Interpreter','latex');
set(h,'rotation',19)
 
 
% load IntSO.txt
% IntSO=IntSO./norm(yexact,inf);
% h=[0.0002; 0.000125; 0.0001; 0.00005; 0.000025];
% figure(4)
% loglog(h,IntSO(1,:),'-sk',h,IntSO(2,:),'-*k',h,IntSO(3,:),'-ok')
% xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
% ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
% legend('EMR','TR','Strang',2);
% axis([2e-5 2.5e-4 3e-7 8e-3]);
% hold on
% plot([4.5e-5 1.8e-4],[0.000015188*10^(0.1) 0.00024*10^(0.1)],'k--')
% h=text(2.9e-5,6.7*10^(-5.9),'Slope=2','FontSize',16,'Interpreter','latex');
% set(h,'rotation',20)

load BlowupIVPFOSp.txt 
h=[0.00005, 0.000025, 0.0000125, 0.000008, 0.000005];
figure(1)
loglog(h,BlowupIVPFOSp(1,:),'-sk',h,BlowupIVPFOSp(2,:),'-*k',h,BlowupIVPFOSp(3,:),'-ok',h,BlowupIVPFOSp(4,:),'-dk')
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('FE','FEA','SpFEA','SpFE',2);
axis([4e-6 6.5e-5 1e-3 2e0]);
hold on
plot([1e-5 4.5e-5],[0.025 0.1125],'k--')
h=text(6e-6,1.5e-2,'Slope=1','FontSize',16,'Interpreter','latex');
set(h,'rotation',17)

% load BlowupIVPSOSp.txt
% h=[0.00005, 0.000025, 0.0000125, 0.000008, 0.000005];
% figure(2)
% loglog(h,BlowupIVPSOSp(1,:),'-sk',h,BlowupIVPSOSp(2,:),'-*k',h,BlowupIVPSOSp(3,:),'-ok',h,BlowupIVPSOSp(4,:),'-dk')
% xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
% ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
% legend('EMR','TR','SoSpFE','Strang',2);
% axis([4e-6 6.5e-5 1e-9 3e-2]);
% hold on
% plot([1e-5 4e-5],[3.125e-6 5e-5],'k--')
% h=text(6e-6,1.1e-6,'Slope=2','FontSize',16,'Interpreter','latex');
% set(h,'rotation',15)

p=4;
lambda=-1;
lamp=lambda*(1-p);
y0=2;            
t0=0;
tf=0.044;
yexact=(((y0^(1-p)+(1/lambda))*exp(lamp*tf))-(1/lambda))^(1/(1-p));

load BlowupIVPFOSp.txt 
BlowupIVPFOSp=BlowupIVPFOSp./yexact;
h=[0.00005, 0.000025, 0.0000125, 0.000008, 0.000005];
figure(3)
loglog(h,BlowupIVPFOSp(1,:),'-sk',h,BlowupIVPFOSp(2,:),'-*k',h,BlowupIVPFOSp(3,:),'-ok',h,BlowupIVPFOSp(4,:),'-dk')
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('FE','FEA','SpFEA','SpFE',2);
axis([4e-6 6.5e-5 1e-4 2e-1]);
hold on
plot([1e-5 4.5e-5],[0.0025 0.01125],'k--')
h=text(6e-6,1.5e-3,'Slope=1','FontSize',16,'Interpreter','latex');
set(h,'rotation',17)
 
 
% load BlowupIVPSOSp.txt
% BlowupIVPSOSp=BlowupIVPSOSp./yexact;
% h=[0.00005, 0.000025, 0.0000125, 0.000008, 0.000005];
% figure(4)
% loglog(h,BlowupIVPSOSp(1,:),'-sk',h,BlowupIVPSOSp(2,:),'-*k',h,BlowupIVPSOSp(3,:),'-ok',h,BlowupIVPSOSp(4,:),'-dk')
% xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
% ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
% legend('EMR','TR','SoSpFE','Strang',2);
% axis([4e-6 6.5e-5 1e-10 3e-3]);
% hold on
% plot([1e-5 4e-5],[3.125e-7 5e-6],'k--')
% h=text(6e-6,1.1e-7,'Slope=2','FontSize',16,'Interpreter','latex');
% set(h,'rotation',15)
clear all;
format long;
%--------------------------------------------------------------------------
%   This program compares second order B-methods with corresponding standard
%   integration methods. In total four methods are used to solve the PDE
%               u_t=u_xx+delta*F(u),  F(u)=(u+alpha)^p.                     (1)
%   The four methods are:
%       EMR    - The explicit midpoint rule. 
%       TR     - The trapezoidal rule.
%       SoSpFE - Second order splitting method constructed with FE.
%       Strang - Second order method constructed from Strang splitting.
%   Each method solves the PDE (1) with 5 different step sizes and 4 
%   different spatially discretizations. This code is parallelized using a 
%   Single Program Multiple Data (SPMD) structure. The explicit methods 
%   have a smaller CPU time, therefore these are carried out serially on 
%   one core while the implicit methods are parallelized on 10 cores.
%
%   NOTE: there needs to be 11 matlab workers in the matlabpool.
%
%--------------------------------------------------------------------------

if matlabpool('size')==0          % if no matlab pool exists.
    matlabpool('open',11);
elseif matlabpool('size')~=11     % if matlab pool does not have 11 cores.
    matlabpool('close');
    matlabpool('open',11);
end

tic                        % start timer.
t0=0;                      % initial time.
tf=0.1150;                 % final time of integration.
n=[25; 51; 75; 101];    % number of spatial discretization nodes.
delta=3;
p=2;
alpha=2;
h=[0.0002; 0.0001; 0.00005; 0.000025; 0.0000125];  % time step sizes.

opt=odeset('AbsTol',1e-12,'RelTol',1e-12);  
parfor i=1:4       % compute reference solution for each n.
    x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].        
    hx=2/(n(i)-1);            % uniform distance between spatial nodes.
    y0=cos(pi*x/2);           % initial condition.
    sol=ode45(@(t,y) fbeck(t,y,n(i),hx,delta,alpha,p),[t0 tf],y0,opt); 
    ytmp=deval(sol,tf);
    yexact{i}=ytmp;
end

spmd(11)   % run single program on 11 cores.  
    [GEEMR, GETR, GE, GESTR]=SoMethods(t0,tf,h,n,p,delta,alpha,yexact);
end

GETR=gather(GETR);                % gather codistributed arrays.
GE=gather(GE);                    
geTR=GETR(:,1:5)+GETR(:,6:10);    % combine results from cores 1-5 and 6-10.
ge=GE(:,1:5)+GE(:,6:10);          
geEMR=GEEMR{11};                  % get results ran on the single core 11.
geStr=GESTR{11};


for i=1:4         % print results to a file.
FileName=strcat('~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE',...
    '/BlowupPDESOSp',num2str(n(i)),'.txt');
fp=fopen(FileName,'w+');
fprintf(fp,'EMR & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geEMR(i,1),geEMR(i,2),geEMR(i,3),geEMR(i,4),geEMR(i,5));
fprintf(fp,'TR & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geTR(i,1),geTR(i,2),geTR(i,3),geTR(i,4),geTR(i,5));
fprintf(fp,'SoSpFE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    ge(i,1),ge(i,2),ge(i,3),ge(i,4),ge(i,5));
fprintf(fp,'Strang & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geStr(i,1),geStr(i,2),geStr(i,3),geStr(i,4),geStr(i,5));
fclose(fp);
end

f1=figure(1);          % create convergence plots.
loglog(h,geEMR(1,:),'-s',h,geTR(1,:),'-*g',h,ge(1,:),'-or',h,geStr(1,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e0]);
print(f1,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE/BlowupPDESO25.eps');

f2=figure(2);
loglog(h,geEMR(2,:),'-s',h,geTR(2,:),'-*g',h,ge(2,:),'-or',h,geStr(2,:),'-dk');
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e-1]);
hold on
plot([3e-5 1.5e-4],[4.5e-5 1.1e-3],'k--')
h=text(1.7e-5,1.6e-5,'Slope=2','FontSize',16,'Interpreter','latex');
set(h,'rotation',19)
print(f2,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE/BlowupPDESO51.eps');

f3=figure(3);
loglog(h,geEMR(3,:),'-s',h,geTR(3,:),'-*g',h,ge(3,:),'-or',h,geStr(3,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e0]);
print(f3,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE/BlowupPDESO75.eps');

f4=figure(4);
loglog(h,geEMR(4,:),'-s',h,geTR(4,:),'-*g',h,ge(4,:),'-or',h,geStr(4,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e0]);
print(f4,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE/BlowupPDESO101.eps');

matlabpool('close');  
toc   % end timer.

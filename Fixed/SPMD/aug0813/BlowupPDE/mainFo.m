clear all;
format long;
%--------------------------------------------------------------------------
%   This program compares first order B-methods with corresponding standard
%   integration methods. In total four methods are used to solve the PDE
%               u_t=u_xx+F(u),  F(u)=(u+alpha)^p.                      (1)
%   The four methods are:
%       FE    - The forward Euler method. 
%       SpFE  - Splitting B-method constructed using FE.
%       FEA   - The backward Euler method, which is the adjoint of FE.
%       SpFEA - Splitting B-method constructed using FEA.
%   Each method solves the IVP (1) with 5 different step sizes and 4 
%   different spatially discretizations. This code is parallelized using a 
%   Single Program Multiple Data (SPMD) structure. The explicit methods 
%   have a smaller CPU time, therefore these are carried out serially on 
%   one core while the implicit methods are parallelized on 10 cores.
%
%   NOTE: there needs to be 11 matlab workers in the matlabpool.
%
%--------------------------------------------------------------------------

if matlabpool('size')==0         % if no matlab pool exists.
    matlabpool('open',11);
elseif matlabpool('size')~=11    % if matlab pool does not have 11 cores.
    matlabpool('close');
    matlabpool('open',11);
end

tic                      % start timer.
t0=0;                    % initial time.
tf=0.1150;               % final time of integration.
n=[25; 51; 75; 101;];  % number of spatial discretization nodes.
delta=3;
p=2;
alpha=2;
h=[0.0001; 0.00005; 0.000025; 0.0000125; 0.000008];  % time step sizes.

opt=odeset('AbsTol',1e-12,'RelTol',1e-12);
parfor i=1:4         % compute reference solution for each n.
    x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
    hx=2/(n(i)-1);            % uniform distance between spatial nodes.
    y0=cos(pi*x/2);           % initial condition.
    sol=ode45(@(t,y) fbeck(t,y,n(i),hx,delta,alpha,p),[t0 tf],y0,opt); 
    ytmp=deval(sol,tf);
    yexact{i}=ytmp;
end

spmd(11)        % run single program on 11 cores. 
   [GEA, GEFEA, GE, GEFE]=FoMethods(t0,tf,h,n,p,delta,alpha,yexact);
end

GEA=gather(GEA);                  % gather codistributed arrays.
GEFEA=gather(GEFEA);
geA=GEA(:,1:5)+GEA(:,6:10);       % combine results from cores 1-5 and 6-10.
geFEA=GEFEA(:,1:5)+GEFEA(:,6:10);
ge=GE{11};                        % get results ran on the single core 11.
geFE=GEFE{11};

for i=1:4             % print results to a file.
FileName=strcat('~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE',...
    '/BlowupPDEFOSp',num2str(n(i)),'.txt');
fp=fopen(FileName,'w+');
fprintf(fp,'FE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geFE(i,1),geFE(i,2),geFE(i,3),geFE(i,4),geFE(i,5));
fprintf(fp,'FEA & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geFEA(i,1),geFEA(i,2),geFEA(i,3),geFEA(i,4),geFEA(i,5));
fprintf(fp,'SpFEA & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geA(i,1),geA(i,2),geA(i,3),geA(i,4),geA(i,5));
fprintf(fp,'SpFE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    ge(i,1),ge(i,2),ge(i,3),ge(i,4),ge(i,5));
fclose(fp);
end

f1=figure(1);         % create convergence plots.
loglog(h,ge(1,:),'-s',h,geFE(1,:),'-*g',h,geA(1,:),'-or',h,geFEA(1,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('SpFE','FE','SpFEA','FEA',2);
axis([5e-6 2e-4 1e-3 1e1]);
print(f1,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE/BlowupPDEFO25.eps');

f2=figure(2);
loglog(h,ge(2,:),'-s',h,geFE(2,:),'-*g',h,geA(2,:),'-or',h,geFEA(2,:),'-dk');
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('FE','FEA','SpFEA','SpFE',2);
axis([5e-6 2e-4 1e-3 1e1]);
hold on
plot([2e-5 8e-5],[0.08 0.32],'k--')
h=text(1e-5,4e-2,'Slope=1','FontSize',16,'Interpreter','latex');
set(h,'rotation',17)
print(f2,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE/BlowupPDEFO51.eps');

f3=figure(3);
loglog(h,ge(3,:),'-s',h,geFE(3,:),'-*g',h,geA(3,:),'-or',h,geFEA(3,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('SpFE','FE','SpFEA','FEA',2);
axis([5e-6 2e-4 1e-3 1e1]);
print(f3,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE/BlowupPDEFO75.eps');

f4=figure(4);
loglog(h,ge(4,:),'-s',h,geFE(4,:),'-*g',h,geA(4,:),'-or',h,geFEA(4,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('SpFE','FE','SpFEA','FEA',2);
axis([5e-6 2e-4 1e-3 1e1]);
print(f4,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/Results/BlowupPDE/BlowupPDEFO101.eps');

matlabpool('close');
toc   % end timer.

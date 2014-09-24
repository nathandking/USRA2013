clear all;
format long;
%--------------------------------------------------------------------------
%   This program compares second order B-methods with corresponding standard
%   integration methods. In total four methods are used to solve the PDE
%               u_t=u_xx+F(u),  F(u)=exp(u).                          (1)
%   The four methods are:
%       EMR    - The explicit midpoint rule. 
%       TR     - The trapezoidal rule.
%       SoSpFE - Second order splitting method constructed with FE.
%       Strang - Second order method constructed from Strang splitting.
%   Each method solves the IVP (1) with 5 different step sizes and 4 
%   different spatially discretizations. This code is parallelized using a 
%   Single Program Multiple Data (SPMD) structure. The explicit methods 
%   have a smaller CPU time, therefore these are carried out serially on 
%   one core while the implicit methods are parallelized on 10 cores.
%
%   NOTE: there needs to be 11 matlab workers in the matlabpool.
%
%--------------------------------------------------------------------------

% open matlab pool of the right size.
if matlabpool('size')==0 
    matlabpool('open',11);
elseif matlabpool('size')~=11
    matlabpool('close');
    matlabpool('open',11);
end

tic    % start timer.
t0=0;                  % initial time.
tf=.19;             % final time of integration.
n=[25; 51; 75; 101];                  % number of spatial nodes.
delta=3;
h=[0.0002; 0.000125; 0.0001; 0.00005; 0.000025];  % time step sizes.

% compute the reference solution.
opt=odeset('AbsTol',1e-12,'RelTol',1e-12);
parfor i=1:4
    x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
    hx=2/(n(i)-1);            % uniform distance between spatial nodes.
    y0=cos(pi*x/2);           % initial condition.
    aij=a(x,x);
    sol=ode45(@(t,y) f(t,y,n(i),hx,delta,aij),[t0 tf],y0,opt); 
    ytmp=deval(sol,tf);
    yexact{i}=ytmp;
end

% run the explicit methods in parallel for the different time steps.
spmd(11)
    [GEEMR, GETR, GESTR]=SoMethods(t0,tf,h,n,delta,yexact);
end
% reconstruct data appropriately.
GETR=gather(GETR);
GESTR=gather(GESTR);
geTR=GETR(:,1:5)+GETR(:,6:10);
geStr=GESTR(:,1:5)+GESTR(:,6:10);
geEMR=GEEMR{11};


% print results to a file.
for i=1:4
FileName=strcat('../../../../../../../tex/USRA2013/IntDiff/Results/IntEqn',...
    '/IntEqnSOSp',num2str(n(i)),'.txt');
fp=fopen(FileName,'w+');
fprintf(fp,'EMR & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geEMR(i,1),geEMR(i,2),geEMR(i,3),geEMR(i,4),geEMR(i,5));
fprintf(fp,'TR & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geTR(i,1),geTR(i,2),geTR(i,3),geTR(i,4),geTR(i,5));
fprintf(fp,'Strang & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',...
    geStr(i,1),geStr(i,2),geStr(i,3),geStr(i,4),geStr(i,5));
fclose(fp);
end

% create convergence plots.
f1=figure(1);
loglog(h,geEMR(1,:),'-s',h,geTR(1,:),'-*g',h,geStr(1,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('EMR','TR','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e-1]);
print(f1,'-deps','../../../../../../../tex/USRA2013/IntDiff/Results/IntEqn/IntEqnSO25.eps');

f2=figure(2);
loglog(h,geEMR(2,:),'-s',h,geTR(2,:),'-*g',h,geStr(2,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('EMR','TR','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e-1]);
print(f2,'-deps','../../../../../../../tex/USRA2013/IntDiff/Results/IntEqn/IntEqnSO51.eps');

f3=figure(3);
loglog(h,geEMR(3,:),'-s',h,geTR(3,:),'-*g',h,geStr(3,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('EMR','TR','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e-1]);
print(f3,'-deps','../../../../../../../tex/USRA2013/IntDiff/Results/IntEqn/IntEqnSO75.eps');

f4=figure(4);
loglog(h,geEMR(4,:),'-s',h,geTR(4,:),'-*g',h,geStr(4,:),'-dk');
xlabel('$\Delta \hspace{0.25cm} t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('EMR','TR','Strang',2);
axis([1e-5 2.5e-4 1e-7 1e-1]);
print(f4,'-deps','../../../../../../../tex/USRA2013/IntDiff/Results/IntEqn/IntEqnSO101.eps');

matlabpool('close');
toc   % end timer.
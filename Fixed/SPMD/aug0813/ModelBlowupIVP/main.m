clear all;
format long;
%--------------------------------------------------------------------------
%   This program compares first order B-methods with corresponding standard
%   integration methods. In total four methods are used to solve the IVP
%               y'=lambda*y+y^p,  y(0)=y0.                          (1)
%   The four methods are:
%       FE    - The forward Euler method. 
%       SpFE  - Splitting B-method constructed using FE.
%       FEA   - The backward Euler method, which is the adjoint of FE.
%       SpFEA - Splitting B-method constructed using FEA.
%   Each method solves the IVP (1) with 5 different step sizes. This code 
%   is therefore parallelized using a Single Program Multiple Data (SPMD)
%   structure. NOTE: there needs to be 5 matlab workers in the matlabpool.
%--------------------------------------------------------------------------

tic   % start timer.
p=4;
lambda=-1;
lamp=lambda*(1-p);
y0=2;            
t0=0;
tf=0.044;
Tb=(1/lamp)*log(1/(lambda*y0^(1-p)+1));  % theoretical blow-up time.
h=[0.00005, 0.000025, 0.0000125, 0.000008, 0.000005];  % time step sizes.
yexact=(((y0^(1-p)+(1/lambda))*exp(lamp*tf))-(1/lambda))^(1/(1-p));

% run the explicit methods in parallel for the different time steps.
spmd(5)
    AA=BlowupODEForwardEuler(t0,tf,y0,h,p,lambda,yexact);
end
% convert output from cell to matrix to allow indexing.
A=cell2mat(AA(:));
ge=A(:,1);
geFE=A(:,2);

% run the implicit methods in parallel for the different time steps.
spmd(5)
    BB=BlowupODEBackwardEuler(t0,tf,y0,h,p,lambda,yexact);
end
% convert output from cell to matrix to allow indexing.
B=cell2mat(BB(:));
geA=B(:,1);
geFEA=B(:,2);

% print results to a file.
FileName='../../../../../../../tex/USRA2013/IntDiff/Results/BlowupIVP/BlowupIVPFOSp.txt';
fp=fopen(FileName,'w+');
fprintf(fp,'FE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geFE(1),...
    geFE(2),geFE(3),geFE(4),geFE(5));
fprintf(fp,'FEA & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geFEA(1),...
    geFEA(2),geFEA(3),geFEA(4),geFEA(5));
fprintf(fp,'SpFEA & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geA(1),...
    geA(2),geA(3),geA(4),geA(5));
fprintf(fp,'SpFE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',ge(1),...
    ge(2),ge(3),ge(4),ge(5));
fclose(fp);

% create convergence plot.
f1=figure(1);
loglog(h,ge,'-s',h,geFE,'-*g',h,geA,'-or',h,geFEA,'-dk');
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('SpFE','FE','SpFEA','FEA',2);
axis([4e-6 6.5e-5 1e-3 2e0]);
hold on
plot([1e-5 4.5e-5],[0.025 0.1125],'k--')
h=text(6e-6,1.5e-2,'Slope=1','FontSize',16,'Interpreter','latex');
set(h,'rotation',17)
print(f1,'-deps','../../../../../../../tex/USRA2013/IntDiff/Results/BlowupIVP/BlowupIVPFO.eps');

toc   % end timer.
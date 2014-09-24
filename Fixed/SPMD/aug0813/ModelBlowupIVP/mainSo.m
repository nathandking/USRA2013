clear all;
format long;
%--------------------------------------------------------------------------
%   This program compares second order B-methods with standard
%   integration methods. In total four methods are used to solve the IVP
%               y'=lambda*y+y^p,  y(0)=y0.                          (1)
%   The four methods are:
%       EMR    - The explicit midpoint rule. 
%       TR     - The trapezoidal rule for differential equations.
%       SoSpFE - The second order splitting B-method constructed from FE.
%       Strang - Strang splitting construction.
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

% run the second order methods in parallel for the different time steps.
spmd(5)
    AA=SoMethods(t0,tf,y0,h,p,lambda,yexact);
end
% convert output from cell to matrix to allow indexing.
A=cell2mat(AA(:));
ge=A(:,1);
geStrang=A(:,2);
geExpMid=A(:,3);
geTrap=A(:,4);

% print results to a file.
fp=fopen('../../../../../../../tex/USRA2013/IntDiff/Results/BlowupIVP/BlowupIVPSOSp.txt','w+');
fprintf(fp,'EMR & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geExpMid(1),geExpMid(2),geExpMid(3),geExpMid(4),geExpMid(5));
fprintf(fp,'TR & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geTrap(1),geTrap(2),geTrap(3),geTrap(4),geTrap(5));
fprintf(fp,'SoSpFE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',ge(1),ge(2),ge(3),ge(4),ge(5));
fprintf(fp,'Strang & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geStrang(1),geStrang(2),geStrang(3),geStrang(4),geStrang(5));
fclose(fp);

% create convergence plot.
f1=figure(1);
loglog(h,geExpMid,'-s',h,geTrap,'-xg',h,ge,'-dk',...
    h,geStrang,'-or');
xlabel('$\log_{\; 10} (\Delta t)$','Interpreter','latex','FontSize',16);
ylabel('$\log_{\; 10} (\mbox{Error })$','Interpreter','latex','FontSize',16);
legend('EMR','TR','SoSpFE','Strang',2)
axis([4e-6 6.5e-5 1e-9 3e-2]);
hold on
plot([1e-5 4e-5],[3.125e-6 5e-5],'k--')
h=text(6e-6,1.1e-6,'Slope=2','FontSize',16,'Interpreter','latex');
set(h,'rotation',15)
print(f1,'-deps','../../../../../../../tex/USRA2013/IntDiff/Results/BlowupIVP/BlowupIVPSO.eps');

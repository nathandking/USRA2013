clear all;
format long;
global p lambda lamp;
p=4;
lambda=-1;
lamp=lambda*(1-p);
y0=2;            % initial condition.
t0=0;
tf=0.044;
Tb=(1/lamp)*log(1/(lambda*y0^(1-p)+1));  % theoretical blow-up time.
ht=[0.00005 0.000025 0.0000125 0.000008 0.000005];  % time step sizes.

[fp fpPlot]=BlowupODEForwardEuler(t0,tf,y0,ht);
 BlowupODEBackwardEuler(fp,fpPlot,t0,tf,y0,ht);


files = dir('*.dat');
for i=1:length(files)
    eval(['load ' files(i).name ' -ascii']);
end

% create convergence plot.
f1=figure(1);
loglog(ht,convplot(1,:),'-s',ht,convplot(2,:),'-*g',...
    ht,convplot(3,:),'-or',ht,convplot(4,:),'-dk');
xlabel('$\Delta t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('SpFE','FE','SpFEA','FEA',2)
axis([2.5e-6 1e-4 1e-3 1e0]);
print(f1,'-deps','~/mun/Dropbox/Haynes-King/tex/IntDiff/BlowupIVPFO.eps');


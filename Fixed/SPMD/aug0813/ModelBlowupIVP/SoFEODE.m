clear all;
format long;
tic
global p lambda;
p=4;
lambda=-1;
lamp=lambda*(1-p);

y0=2;      % initial condition.
ybeck(1)=y0;
y(1)=y0;
z(1)=y0;
Tb=(1/lamp)*log(1/(lambda*y0^(1-p)+1));  % theoretical blow-up time.
  
% Initial conditions.
t(1)=0;
ht=[0.00005,0.000025,0.0000125,0.000008,0.000005];
for j=1:5
t=0:ht(j):0.044;
    
   for i=1:length(t)
    % Integrate one step.
    %Explicit Midpoint.
    zhat(i+1)=z(i)+0.5*ht(j)*(lambda*z(i)+z(i)^p);
    z(i+1)=z(i)+ht(j)*(lambda*zhat(i+1)+zhat(i+1)^p);
    
    %Beck.
    ybeck(i+1)=((1-p)*ht(j)/2+(((1+lambda*ht(j)/2)/(1-lambda*ht(j)/2))^(1-p))*((1-p)*ht(j)/2+ybeck(i)^(1-p)))^(1/(1-p));
    
    %Strang Splitting.
    y1=diffus(y(i),ht(j)/2);
    y2=blowup(y1,ht(j));
    y3=diffus(y2,ht(j)/2);
    y(i+1)=y3;
   
    yexact(i)=(((y0^(1-p)+(1/lambda))*exp(lamp*t(i)))-(1/lambda))^(1/(1-p));
    geplot(i)=abs(y(i)-yexact(i)');
   end
     ge(j)=abs(ybeck(i)-yexact(i)');
     geStrang(j)=abs(y(i)-yexact(i)');
     geExpMid(j)=abs(z(i)-yexact(i)')
end
ge'
geStrang'
geExpMid'
% Printing results
fprintf('Order test:\n')
fprintf('SoSpFE       Strang       Exp Mid\n')
fprintf('%10.8f %10.8f %10.8f\n',ge(1)/ge(2),geStrang(1)/geStrang(2),geExpMid(1)/geExpMid(2));

fp=fopen('~/Dropbox/Haynes-King/tex/IntDiff/BlowupIVPSOSp.txt','w+');
fprintf(fp,'EMR & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geExpMid(1),geExpMid(2),geExpMid(3),geExpMid(4),geExpMid(5));
fprintf(fp,'SoSpFE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',ge(1),ge(2),ge(3),ge(4),ge(5));
fprintf(fp,'Strang & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geStrang(1),geStrang(2),geStrang(3),geStrang(4),geStrang(5));
fclose(fp);

% create convergence plot.
f1=figure(1);
loglog(ht,geExpMid,'-s',ht,ge,'-dk',...
    ht,geStrang,'-or');
xlabel('$\Delta t$','Interpreter','latex','FontSize',16);
ylabel('Error','FontSize',16);
legend('EMR','SoSpFE','Strang',2)
axis([2.5e-6 1e-4 1e-9 1e-2]);
print(f1,'-deps','~/Dropbox/Haynes-King/tex/IntDiff/BlowupIVPSO.eps');

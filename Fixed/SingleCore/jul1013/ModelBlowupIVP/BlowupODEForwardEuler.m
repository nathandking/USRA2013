function [fp,fpPlot]=BlowupODEForwardEuler(t0,tf,y0,ht)
global p lambda lamp;
 tic
% Initial conditions.
yexact=(((y0^(1-p)+(1/lambda))*exp(lamp*tf))-(1/lambda))^(1/(1-p)); 

parfor j=1:5
t=t0:ht(j):tf;
z=y0;
y=y0;
y1=zeros(1,1)
    for i=1:length(t)-1
    % Integrate one step.
%     if i<=2
%     size(y)
%     i
%     size(z)
%     end
    z=z+ht(j)*(lambda*z+z^p);
    k1=y+ht(j)*lambda*y;
    y1=y+(ht(j)/2)*(lambda*y+lambda*k1);
    
    y=((1-p)*ht(j)+y1.^(1-p)).^(1./(1-p))
    %disp(y)
    %size(y)
%     if i==1
%     size(y)
%     y1=diffus(y,ht(j));
%     size(y1)
%     
%     y=blowup(y,ht(j));
%     size(y)
%     end
    end
% size(y)
% size(z)
ge(j)=abs(y-yexact);
geFE(j)=abs(z-yexact);
end
toc
ge'
% Printing results
fprintf('Order test:\n')
fprintf('FE           B-Method Euler\n')
fprintf('%10.8f      %10.8f\n',geFE(1)/geFE(2),ge(1)/ge(2))


FileName=strcat('~/mun/Dropbox/Haynes-King/tex/IntDiff/BlowupIVPFOSp','.txt');
fpPlot=fopen('convplot.dat','w+');
fprintf(fpPlot,'%10.8e %10.8e %10.8e %10.8e %10.8e\n',ge(1),ge(2),ge(3),ge(4),ge(5));
fprintf(fpPlot,'%10.8e %10.8e %10.8e %10.8e %10.8e\n',geFE(1),geFE(2),geFE(3),geFE(4),geFE(5));
fp=fopen(FileName,'w+');
fprintf(fp,'SpFE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',ge(1),ge(2),ge(3),ge(4),ge(5));
fprintf(fp,'FE & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geFE(1),geFE(2),geFE(3),geFE(4),geFE(5));


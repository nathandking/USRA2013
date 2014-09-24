function []=BlowupODEBackwardEuler(fp,fpPlot,t0,tf,y0,ht)
global p lambda lamp;

parfor j=1:5
[getmp, geFEAtmp]=StupidCode(ht(j),t0,tf,y0);
ge(j)=getmp;
geFEA(j)=geFEAtmp;
end
ge'
% Printing results
fprintf('Order test:\n')
fprintf('FEA           B-Method Backward Euler\n')
fprintf('%10.8f      %10.8f\n',geFEA(1)/geFEA(2),ge(1)/ge(2))

fprintf(fpPlot,'%10.8e %10.8e %10.8e %10.8e %10.8e\n',ge(1),ge(2),ge(3),ge(4),ge(5));
fprintf(fpPlot,'%10.8e %10.8e %10.8e %10.8e %10.8e\n',geFEA(1),geFEA(2),geFEA(3),geFEA(4),geFEA(5));
fprintf(fp,'SpFEA & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',ge(1),ge(2),ge(3),ge(4),ge(5));
fprintf(fp,'FEA & %6.3e & %6.3e & %6.3e & %6.3e & %6.3e\\\\\n',geFEA(1),geFEA(2),geFEA(3),geFEA(4),geFEA(5));
fclose(fp);
fclose(fpPlot);
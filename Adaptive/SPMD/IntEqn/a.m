function [axy]=a(x,y)

for i=1:length(x)
    for j=1:length(y)
        if (x(i)>y(j))
            axy(i,j)=exp(-(x(i)-y(j)));
        elseif (x(i)<y(j))
            axy(i,j)=exp(-(1+x(i)-y(j)));
        elseif (x(i)==y(j))
            axy(i,j)=0.5*(1+exp(-1));
        end
    end
end
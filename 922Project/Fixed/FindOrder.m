   
err=[0.003331980135727; 0.000841838454562; 0.000220013510329; 0.000066680911004];
hx=[1e-04; 5e-05; 2.5e-05; 1.25e-05];
%hx=[1e-04; 5e-05; 2.5e-05; 1.25e-05; 8e-06];
%hx=[2e-04 ;  1.25e-04 ;  1e-04 ;  5e-05 ;  2.5e-05];
%hx=[ 5e-05 ;  2.5e-05 ;  1.25e-05 ;  8e-06 ;  5e-06];

p=zeros(3,1);
for i=1:3
    p(i)=(log(err(i))-log(err(i+1)))/(log(hx(i))-log(hx(i+1)));
end

pavg=sum(p)/3;
fprintf('Average order = %8.4f\n',pavg)
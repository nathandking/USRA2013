err=[3.312e-04; 1.641e-06; 3.403e-07; 7.883e-08; 1.903e-08];
hx=[2e-04; 1e-04; 5e-05; 2.5e-05; 1.25e-05];
%hx=[1e-04; 5e-05; 2.5e-05; 1.25e-05; 8e-06];
%hx=[2e-04 ;  1.25e-04 ;  1e-04 ;  5e-05 ;  2.5e-05];
%hx=[ 5e-05 ;  2.5e-05 ;  1.25e-05 ;  8e-06 ;  5e-06];

p=zeros(4,1);
for i=1:4
    p(i)=(log(err(i))-log(err(i+1)))/(log(hx(i))-log(hx(i+1)));
end

pavg=sum(p)/4;
fprintf('Average order = %8.4f\n',pavg)

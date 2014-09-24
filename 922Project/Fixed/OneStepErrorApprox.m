t(1)=0;                % initial time.
tf=0.493;               % final time of integration. (Tb=4.944521e-01)
n=51;                  % number of spatial nodes.
x=linspace(-1,1,n);    % define spatial nodes on the inteval [-1,1].          
hx=2/(n-1);            % uniform distance between spatial nodes.
alpha=1;
beta=4;

ySpBE(:,1)=beta*cos(0.5*pi*x);    % initial condition.

options=optimset('TolFun',1e-10);
opt=odeset('AbsTol',1e-12,'RelTol',1e-12);
[T yexact]=ode15s(@f,[t(1) tf],ySpBE(:,1),opt);  % reference solution.

Un=yexact(end,:);
diffUn=zeros(1,n);
diffUn(2:n-1)=(Un(1:n-2)-2*Un(2:n-1)+Un(3:n))/hx^2;
diffdiffUn=zeros(1,n);
diffdiffUn(2:n-1)=(diffUn(1:n-2)-2*diffUn(2:n-1)+diffUn(3:n))/hx^2;
diffUnalpha=zeros(1,n);
diffUnalpha(2:n-1)=((Un(1:n-2)).^(1+alpha)-2*(Un(2:n-1)).^(1+alpha)+(Un(3:n)).^(1+alpha))/hx^2;

k=0.0000125;
SpBE=norm(0.5*k^2*((Un).^(1+alpha).*diffUn- diffUnalpha-diffdiffUn),inf)

BE=norm(-0.5*k^2*((Un).^(1+alpha).*diffUn+diffUnalpha +diffdiffUn + (Un).^(2+2*alpha)),inf) 
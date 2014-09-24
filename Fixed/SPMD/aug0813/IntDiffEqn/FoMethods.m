function [GEA, GEFEA, ge, geFE]=FoMethods(t0,tf,h,n,delta,yexact)

% allocate memory for global error arrays.
ge=zeros(length(n),length(h));
geFE=zeros(length(n),length(h));
geA=zeros(length(n),1);
geFEA=zeros(length(n),1);

if (labindex <= 5)
    for i=4:-2:2            % run implicit methods in parallel.
        ht=h(labindex);
        t=t0:ht:tf;            % define vector of time steps.
        x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
        hx=2/(n(i)-1);            % uniform distance between spatial nodes.
        y0=cos(pi*x/2);        % initial condition.
        y=zeros(n(i),length(t));
        yFEA=zeros(n(i),length(t));
        y(:,1)=y0;
        yFEA(:,1)=y0;
        aij=a(x,x);
        options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,...
            'Display','off');  
        A=(-2*diag(ones(n(i)-2,1))+diag(ones(n(i)-3,1),1)+...
            diag(ones(n(i)-3,1),-1))/hx^2;
        B=(eye(n(i)-2)-ht*A);
    
        for k=1:length(t)
            % integrate one step.
            % backward Euler.
            yFEA(:,k+1)=fsolve(@(x) FEA(x,yFEA(:,k),ht,hx,n(i),delta,aij)...
                ,yFEA(:,k),options);
        
            % B-method.
            y1=B\y(2:n(i)-1,k);
            
            y2=zeros(n(i),1);
            y2(2:n(i)-1)=-log(exp(-y1)-delta*ht);
            
            y(:,k+1)=intdiff(y2(:),ht,hx,n(i),aij);
        end
    geA(i)=norm(y(:,k)-yexact{i},inf);
    geFEA(i)=norm(yFEA(:,k)-yexact{i},inf);
    end
elseif (labindex > 5) && (labindex <= 10)
    for i=3:-2:1            % run implicit methods in parallel
        ht=h(labindex-5);
        t=t0:ht:tf;            % define vector of time steps.
        x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
        hx=2/(n(i)-1);            % uniform distance between spatial nodes.
        y0=cos(pi*x/2);        % initial condition.
        y=zeros(n(i),length(t));
        yFEA=zeros(n(i),length(t));
        y(:,1)=y0;
        yFEA(:,1)=y0;
        aij=a(x,x);
        options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,...
            'Display','off');
        A=(-2*diag(ones(n(i)-2,1))+diag(ones(n(i)-3,1),1)+...
            diag(ones(n(i)-3,1),-1))/hx^2;
        B=(eye(n(i)-2)-ht*A);

        for k=1:length(t)
            % integrate one step.
            % backward Euler.
            yFEA(:,k+1)=fsolve(@(x) FEA(x,yFEA(:,k),ht,hx,n(i),delta,aij)...
                ,yFEA(:,k),options);
        
            % B-method.
            y1=B\y(2:n(i)-1,k);
            
            y2=zeros(n(i),1);
            y2(2:n(i)-1)=-log(exp(-y1)-delta*ht);
            
            y(:,k+1)=intdiff(y2(:),ht,hx,n(i),aij);
        end
    geA(i)=norm(y(:,k)-yexact{i},inf);
    geFEA(i)=norm(yFEA(:,k)-yexact{i},inf);
    end
elseif (labindex == 11)
    for i=4:-1:1               % run explicit first order methods in serial.
        for j=5:-1:1
            t=t0:h(j):tf;            % define vector of time steps.
            x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
            hx=2/(n(i)-1);            % uniform distance between spatial nodes.
            y0=cos(pi*x/2);        % initial condition.
            z=zeros(n(i),length(t));
            yFE=zeros(n(i),length(t));
            z(:,1)=y0;
            yFE(:,1)=y0;
            aij=a(x,x);
            for k=1:length(t)
                % integrate one step.
                % forward Euler.
                yFE(:,k+1)=yFE(:,k)+h(j)*f(k*h(j),yFE(:,k),n(i),...
                    hx,delta,aij);
        
                % B-method.
                z1=zeros(n(i),1);
                z1(2:n(i)-1)=z(2:n(i)-1,k)+h(j)*(z(1:n(i)-2,k)-2*z(2:n(i)-1,k)...
                    +z(3:n(i),k))/hx^2;
        
                z2=zeros(n(i),1);
                z2(2:n(i)-1)=-log(exp(-z1(2:n(i)-1))-delta*h(j));
        
                z(:,k+1)=intdiff(z2(:),h(j),hx,n(i),aij);
            end
            % calculate absolute errors.
            ge(i,j)=norm(z(:,k)-yexact{i},inf);
            geFE(i,j)=norm(yFE(:,k)-yexact{i},inf);
        end
    end
end
codis1=codistributor1d(2,[1 1 1 1 1 1 1 1 1 1 1],[4 11]);
GEA=codistributed.build(geA,codis1);
GEFEA=codistributed.build(geFEA,codis1);

function [GEA, GEFEA, ge, geFE]=FoMethods(t0,tf,h,n,p,delta,alpha,yexact)
%--------------------------------------------------------------------------
%    This function contains the single program to be ran on multiple cores.
%--------------------------------------------------------------------------

ge=zeros(length(n),length(h));   % allocate memory for global error arrays.
geFE=zeros(length(n),length(h));
geA=zeros(length(n),1);
geFEA=zeros(length(n),1);

if (labindex <= 5)
    for i=4:-2:2            % run implicit methods in parallel.
        ht=h(labindex);
        t=t0:ht:tf;            % define vector of time steps.
        x=linspace(-1,1,n(i)); % define spatial nodes on the inteval [-1,1].          
        hx=2/(n(i)-1);         % uniform distance between spatial nodes.
        y0=cos(pi*x/2);        % initial condition.
        y=zeros(n(i),length(t));
        yFEA=zeros(n(i),length(t));
        y(:,1)=y0;
        yFEA(:,1)=y0;
        options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,...
            'Display','off'); 
        A=(-2*diag(ones(n(i)-2,1))+diag(ones(n(i)-3,1),1)+...
            diag(ones(n(i)-3,1),-1))/hx^2;
        B=(eye(n(i)-2)-ht*A);
        for k=1:length(t)         % integrate until tf.
            yFEA(:,k+1)=fsolve(@(x) PsiFEA(x,yFEA(:,k),...    % FEA.
                ht,n(i),hx,delta,alpha,p),yFEA(:,k),options);
        
            y1=B\y(2:n(i)-1,k);      % SpFEA.
            y(2:n(i)-1,k+1)=((y1+alpha).^(1-p)+(1-p)*delta*ht).^...
                (1/(1-p))-alpha;
        end
    % calculate absolute errors.   
    geA(i)=norm(y(:,k)-yexact{i},inf);
    geFEA(i)=norm(yFEA(:,k)-yexact{i},inf);
    end
elseif (labindex > 5) && (labindex <= 10)
    for i=3:-2:1            % run implicit methods in parallel
        ht=h(labindex-5);
        t=t0:ht:tf;            % define vector of time steps.
        x=linspace(-1,1,n(i)); % define spatial nodes on the inteval [-1,1].          
        hx=2/(n(i)-1);         % uniform distance between spatial nodes.
        y0=cos(pi*x/2);        % initial condition.
        y=zeros(n(i),length(t));
        yFEA=zeros(n(i),length(t));
        y(:,1)=y0;
        yFEA(:,1)=y0;
        options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,...
            'Display','off'); 
        A=(-2*diag(ones(n(i)-2,1))+diag(ones(n(i)-3,1),1)+...
            diag(ones(n(i)-3,1),-1))/hx^2;
        B=(eye(n(i)-2)-ht*A);
        for k=1:length(t)         % integrate until tf.
            yFEA(:,k+1)=fsolve(@(x) PsiFEA(x,yFEA(:,k),...    % FEA.
                ht,n(i),hx,delta,alpha,p),yFEA(:,k),options);
        
            y1=B\y(2:n(i)-1,k);      % SpFEA.
            y(2:n(i)-1,k+1)=((y1+alpha).^(1-p)+(1-p)*delta*ht).^...
                (1/(1-p))-alpha;
        end
    % calculate absolute errors.   
    geA(i)=norm(y(:,k)-yexact{i},inf);
    geFEA(i)=norm(yFEA(:,k)-yexact{i},inf);
    end
elseif (labindex == 11)
    for i=4:-1:1               % run explicit first order methods in serial.
        for j=5:-1:1
            t=t0:h(j):tf;            % define vector of time steps.
            x=linspace(-1,1,n(i));   % define spatial nodes on the inteval [-1,1].          
            hx=2/(n(i)-1);           % uniform distance between spatial nodes.
            y0=cos(pi*x/2);          % initial condition.
            z=zeros(n(i),length(t));
            yFE=zeros(n(i),length(t));
            z(:,1)=y0;
            yFE(:,1)=y0;
            for k=1:length(t)          % integrate until tf.
                 yFE(:,k+1)=yFE(:,k)+h(j)*fbeck(k*h(j),...   % FE.
                     yFE(:,k),n(i),hx,delta,alpha,p);
        
                z1=zeros(n(i),1);       % SpFE.
                z1(2:n(i)-1)=((z(2:n(i)-1,k)+alpha*ones(1,n(i)-2)').^...
                    (1-p)+(1-p)*delta*h(j)).^(1/(1-p))-alpha;
                z(2:n(i)-1,k+1)=z1(2:n(i)-1)+h(j)*(z1(1:n(i)-2)-2*...
                    z1(2:n(i)-1)+z1(3:n(i)))/(hx^2);
            end
            % calculate absolute errors.
            ge(i,j)=norm(z(:,k)-yexact{i},inf);
            geFE(i,j)=norm(yFE(:,k)-yexact{i},inf);
        end
    end
end

codis1=codistributor1d(2,[1 1 1 1 1 1 1 1 1 1 1],[4 11]);
GEA=codistributed.build(geA,codis1);      % construct codistributed arrays.
GEFEA=codistributed.build(geFEA,codis1);

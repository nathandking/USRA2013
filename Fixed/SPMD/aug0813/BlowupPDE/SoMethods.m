function [geEMR, GETR, GE, geStr]=SoMethods(t0,tf,h,n,p,delta,alpha,yexact)
%----------------------------------------------------------------------------
%    This function contains the single program to be ran on multiple cores.
%----------------------------------------------------------------------------

ge=zeros(length(n),1);     % allocate memory for global error arrays.
geStr=zeros(length(n),length(h));
geEMR=zeros(length(n),length(h));
geTR=zeros(length(n),1);

if (labindex <= 5)
    for i=4:-2:2            % run implicit methods in parallel.
        ht=h(labindex);
        t=t0:ht:tf;            % define vector of time steps.
        x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
        hx=2/(n(i)-1);            % uniform distance between spatial nodes.
        y0=cos(pi*x/2);        % initial condition.
        y=zeros(n(i),length(t));
        yTR=zeros(n(i),length(t));
        y(:,1)=y0;
        yTR(:,1)=y0;
        options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,'Display','off');
        for k=1:length(t)         % integrate until tf.
            yTR(:,k+1)=fsolve(@(x) trap(x,yTR(:,k),ht,...     % TR.
                n(i),hx,delta,alpha,p),yTR(:,k),options);
        
            v=fsolve(@(x) SoPsi(x,y(:,k),ht,n(i),...        % SpSoFE.
                hx,delta,alpha,p),y(:,k),options);
            y(2:n(i)-1,k+1)=((v(2:n(i)-1)+(0.5*ht/hx^2)*(v(1:n(i)-2)-...
                2*v(2:n(i)-1)+v(3:n(i)))+alpha).^(1-p)+0.5*(1-p)*delta...
                *ht).^(1/(1-p))-alpha;
        end
    % calculate absolute errors.
    ge(i)=norm(y(:,k)-yexact{i},inf);
    geTR(i)=norm(yTR(:,k)-yexact{i},inf);
    end
elseif (labindex > 5) && (labindex <= 10)
    for i=3:-2:1            % run implicit methods in parallel
        ht=h(labindex-5);
        t=t0:ht:tf;            % define vector of time steps.
        x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
        hx=2/(n(i)-1);            % uniform distance between spatial nodes.
        y0=cos(pi*x/2);        % initial condition.
        y=zeros(n(i),length(t));
        yTR=zeros(n(i),length(t));
        y(:,1)=y0;
        yTR(:,1)=y0;
        options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,'Display','off');
        for k=1:length(t)       % integrate until tf.
            yTR(:,k+1)=fsolve(@(x) trap(x,yTR(:,k),ht,n(i),...    % TR.
                hx,delta,alpha,p),yTR(:,k),options);
        
            v=fsolve(@(x) SoPsi(x,y(:,k),ht,n(i),hx,...       % SpSoFE.
                delta,alpha,p),y(:,k),options);
            y(2:n(i)-1,k+1)=((v(2:n(i)-1)+(0.5*ht/hx^2)*(v(1:n(i)-2)-...
                2*v(2:n(i)-1)+v(3:n(i)))+alpha).^(1-p)+0.5*(1-p)*...
                delta*ht).^(1/(1-p))-alpha;
        end
    % calculate absolute errors.
    ge(i)=norm(y(:,k)-yexact{i},inf);
    geTR(i)=norm(yTR(:,k)-yexact{i},inf);
    end
elseif (labindex == 11)
    for i=4:-1:1               % run explicit first order methods in serial.
        for j=5:-1:1
            t=t0:h(j):tf;            % define vector of time steps.
            x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
            hx=2/(n(i)-1);            % uniform distance between spatial nodes.
            y0=cos(pi*x/2);        % initial condition.
            yEMR=zeros(n(i),length(t));
            yStr=zeros(n(i),length(t));
            yEMR(:,1)=y0;
            yStr(:,1)=y0;
            for k=1:length(t)   % integrate until tf.
                z1=yEMR(:,k)+0.5*h(j)*fbeck(k*h(j),...   %EMR.
                    yEMR(:,k),n(i),hx,delta,alpha,p);
                z1(1)=0;
                z1(n(i))=0;
                yEMR(:,k+1)=yEMR(:,k)+h(j)*fbeck(k*h(j),z1,n(i),hx,delta,alpha,p);
                yEMR(1,k+1)=0;
                yEMR(n(i),k+1)=0;
       
                y1=zeros(n(i),1);        % Strang.     
                y1(2:n(i)-1)=((yStr(2:n(i)-1,k)+alpha).^(1-p)+0.5*(1-p)*...
                    delta*h(j)).^(1/(1-p))-alpha;
                y2=ExpMid(y1,h(j),hx,n(i));
                yStr(2:n(i)-1,k+1)=((y2(2:n(i)-1)+alpha).^(1-p)+0.5*(1-p)*...
                    delta*h(j)).^(1/(1-p))-alpha;
            end
            % calculate absolute errors.
            geEMR(i,j)=norm(yEMR(:,k)-yexact{i},inf);
            geStr(i,j)=norm(yStr(:,k)-yexact{i},inf);
        end
    end
end

codis1=codistributor1d(2,[1 1 1 1 1 1 1 1 1 1 1],[4 11]);
GETR=codistributed.build(geTR,codis1);   % construct codistributed arrays.
GE=codistributed.build(ge,codis1);

function [geEMR, GETR, GE, geStr]=SoMethods(t0,tf,h,n,delta,yexact)

% allocate memory for global error arrays.
ge=zeros(length(n),1);
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
    
        for k=1:length(t)
            % integrate one step.
            % trapezoidal rule.
            yTR(:,k+1)=fsolve(@(x) trap(x,yTR(:,k),ht,n(i),hx,delta)...
                ,yTR(:,k),options);
        
            % B-method.
            v1=zeros(n(i),1);
            v1(2:n(i)-1)=-log(exp(-y(2:n(i)-1,k))-0.5*delta*ht);
            
            v2=zeros(n(i),1);
            v2(2:n(i)-1)=v1(2:n(i)-1)+0.5*ht*(v1(1:n(i)-2)-2*v1(2:n(i)-1)+v1(3:n(i)))/hx^2;
            
            v3=fsolve(@(x) ImpEul(x,v2,ht,n(i),hx),v2,...
                options);
        
            y(2:n(i)-1,k+1)=-log(exp(-v3(2:n(i)-1))-0.5*delta*ht);
        end
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
    
        for k=1:length(t)
            % integrate one step.
            % trapezoidal rule.
            yTR(:,k+1)=fsolve(@(x) trap(x,yTR(:,k),ht,n(i),hx,delta)...
                ,yTR(:,k),options);
        
            % B-method.
            v1=zeros(n(i),1);
            v1(2:n(i)-1)=-log(exp(-y(2:n(i)-1,k))-0.5*delta*ht);
            
            v2=zeros(n(i),1);
            v2(2:n(i)-1)=v1(2:n(i)-1)+0.5*ht*(v1(1:n(i)-2)-2*v1(2:n(i)-1)+v1(3:n(i)))/hx^2;
            
            v3=fsolve(@(x) ImpEul(x,v2,ht,n(i),hx),v2,...
                options);
        
            y(2:n(i)-1,k+1)=-log(exp(-v3(2:n(i)-1))-0.5*delta*ht);
        end
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
    
            for k=1:length(t)
                % integrate one step.
                % explicit midpoint rule.
                z1=yEMR(:,k)+0.5*h(j)*fbeck(k*h(j),yEMR(:,k),n(i),hx,delta);
                z1(1)=0;
                z1(n(i))=0;
        
                yEMR(:,k+1)=yEMR(:,k)+h(j)*fbeck(k*h(j),z1,n(i),hx,delta);
                yEMR(1,k+1)=0;
                yEMR(n(i),k+1)=0;
        
                % Strang.
                y1=zeros(n(i),1);
                y1(2:n(i)-1)=-log(exp(-yStr(2:n(i)-1,k))-0.5*delta*h(j));
        
                y2=ExpMid(y1,h(j),hx,n(i));
        
                yStr(2:n(i)-1,k+1)=-log(exp(-y2(2:n(i)-1))-0.5*delta*h(j));
            end
            % calculate absolute errors.
            geEMR(i,j)=norm(yEMR(:,k)-yexact{i},inf);
            geStr(i,j)=norm(yStr(:,k)-yexact{i},inf);
        end
    end
end
codis1=codistributor1d(2,[1 1 1 1 1 1 1 1 1 1 1],[4 11]);
GETR=codistributed.build(geTR,codis1);
GE=codistributed.build(ge,codis1);
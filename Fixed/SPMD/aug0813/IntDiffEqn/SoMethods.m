function [geEMR, GETR, GESTR]=SoMethods(t0,tf,h,n,delta,yexact)

% allocate memory for global error arrays.
geStr=zeros(length(n),1);
geEMR=zeros(length(n),length(h));
geTR=zeros(length(n),1);

if (labindex <= 5)
    for i=4:-2:2            % run implicit methods in parallel.
        ht=h(labindex);
        t=t0:ht:tf;            % define vector of time steps.
        x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
        hx=2/(n(i)-1);            % uniform distance between spatial nodes.
        y0=cos(pi*x/2);        % initial condition.
        yTR=zeros(n(i),length(t));
        yStr=zeros(n(i),length(t));
        yTR(:,1)=y0;
        yStr(:,1)=y0;
        aij=a(x,x);
        options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,'Display','off');
    
        for k=1:length(t)
            % integrate one step.
            % trapezoidal rule.
            yTR(:,k+1)=fsolve(@(x) trap(x,yTR(:,k),ht,n(i),hx,delta,aij)...
                ,yTR(:,k),options);  
            
            % Strang.
            y15=zeros(n(i),1);
            y15(2:n(i)-1)=-log(exp(-yStr(2:n(i)-1,k))-delta*ht/2);
            y25=ExpSecOrd(y15,ht/2,n(i),hx);
            y35=intdiff(y25,ht,hx,n(i),aij);
            y45=ExpSecOrd(y35,ht/2,n(i),hx);
            yStr(2:n(i)-1,k+1)=-log(exp(-y45(2:n(i)-1))-delta*ht/2);
        end
    geTR(i)=norm(yTR(:,k)-yexact{i},inf);
    geStr(i)=norm(yStr(:,k)-yexact{i},inf);
    end
elseif (labindex > 5) && (labindex <= 10)
    for i=3:-2:1            % run implicit methods in parallel
        ht=h(labindex-5);
        t=t0:ht:tf;            % define vector of time steps.
        x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
        hx=2/(n(i)-1);            % uniform distance between spatial nodes.
        y0=cos(pi*x/2);        % initial condition.
        yTR=zeros(n(i),length(t));
        yStr=zeros(n(i),length(t));
        yTR(:,1)=y0;
        yStr(:,1)=y0;
        aij=a(x,x);
        options=optimoptions('fsolve','TolFun',1e-14,'Tolx',1e-14,'Display','off');
    
        for k=1:length(t)
            % integrate one step.
            % trapezoidal rule.
            yTR(:,k+1)=fsolve(@(x) trap(x,yTR(:,k),ht,n(i),hx,delta,aij)...
                ,yTR(:,k),options);
            
            % Strang.
            y15=zeros(n(i),1);
            y15(2:n(i)-1)=-log(exp(-yStr(2:n(i)-1,k))-delta*ht/2);
            y25=ExpSecOrd(y15,ht/2,n(i),hx);
            y35=intdiff(y25,ht,hx,n(i),aij);
            y45=ExpSecOrd(y35,ht/2,n(i),hx);
            yStr(2:n(i)-1,k+1)=-log(exp(-y45(2:n(i)-1))-delta*ht/2);
        end
    geTR(i)=norm(yTR(:,k)-yexact{i},inf);
    geStr(i)=norm(yStr(:,k)-yexact{i},inf);
    end
elseif (labindex == 11)
    for i=4:-1:1               % run explicit first order methods in serial.
        for j=5:-1:1
            t=t0:h(j):tf;            % define vector of time steps.
            x=linspace(-1,1,n(i));    % define spatial nodes on the inteval [-1,1].          
            hx=2/(n(i)-1);            % uniform distance between spatial nodes.
            y0=cos(pi*x/2);        % initial condition.
            yEMR=zeros(n(i),length(t));
            yEMR(:,1)=y0;
            aij=a(x,x);
            for k=1:length(t)
                % integrate one step.
                % explicit midpoint rule.               
                z1=yEMR(:,k)+0.5*h(j)*f(k*h(j),yEMR(:,k),n(i),hx,delta,aij);
                z1(1)=0;
                z1(n(i))=0;
        
                yEMR(:,k+1)=yEMR(:,k)+h(j)*f(k*h(j),z1,n(i),hx,delta,aij);
                yEMR(1,k+1)=0;
                yEMR(n(i),k+1)=0;
            end
            % calculate absolute errors.
            geEMR(i,j)=norm(yEMR(:,k)-yexact{i},inf);
        end
    end
end
codis1=codistributor1d(2,[1 1 1 1 1 1 1 1 1 1 1],[4 11]);
GETR=codistributed.build(geTR,codis1);
GESTR=codistributed.build(geStr,codis1);
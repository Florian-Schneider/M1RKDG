function Y=InitialCondition(x,parameter)

global PVal


Weights = parameter.Domain.Quadrature.Weights;

dx = parameter.Domain.dx;
dy = parameter.Domain.dy;

QuadraturePointsx = parameter.Domain.Quadrature.PhysicalPointsx;
QuadraturePointsy = parameter.Domain.Quadrature.PhysicalPointsy;

p1 = @(x,y) 1;
p2 = @(x,y) x;
p3 = @(x,y) y;
p4 = @(x,y) p2(x,y).*p3(x,y);
p5 = @(x,y) p2(x,y).^2-1/3;
p6 = @(x,y) p3(x,y).^2-1/3;

p2ValLeft = p2(parameter.Domain.Quadrature.LeftEdgex,parameter.Domain.Quadrature.LeftEdgey);
p2ValRight = p2(parameter.Domain.Quadrature.RightEdgex,parameter.Domain.Quadrature.RightEdgey);
p2ValUpper = p2(parameter.Domain.Quadrature.UpperEdgex,parameter.Domain.Quadrature.UpperEdgey);
p2ValLower = p2(parameter.Domain.Quadrature.LowerEdgex,parameter.Domain.Quadrature.LowerEdgey);
p2Val = p2(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p2Gradx = 2/dx+0*p2Val;
p2Grady = 0*p2Val;

p3ValLeft = p3(parameter.Domain.Quadrature.LeftEdgex,parameter.Domain.Quadrature.LeftEdgey);
p3ValRight = p3(parameter.Domain.Quadrature.RightEdgex,parameter.Domain.Quadrature.RightEdgey);
p3ValUpper = p3(parameter.Domain.Quadrature.UpperEdgex,parameter.Domain.Quadrature.UpperEdgey);
p3ValLower = p3(parameter.Domain.Quadrature.LowerEdgex,parameter.Domain.Quadrature.LowerEdgey);
p3Val = p3(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p3Gradx = 0*p3Val;
p3Grady = 2/dy+0*p3Val;

p4ValLeft = p4(parameter.Domain.Quadrature.LeftEdgex,parameter.Domain.Quadrature.LeftEdgey);
p4ValRight = p4(parameter.Domain.Quadrature.RightEdgex,parameter.Domain.Quadrature.RightEdgey);
p4ValUpper = p4(parameter.Domain.Quadrature.UpperEdgex,parameter.Domain.Quadrature.UpperEdgey);
p4ValLower = p4(parameter.Domain.Quadrature.LowerEdgex,parameter.Domain.Quadrature.LowerEdgey);
p4Val = p4(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p4Gradx = 2/dx*p3Val;
p4Grady = 2/dy*p2Val;

p5ValLeft = p5(parameter.Domain.Quadrature.LeftEdgex,parameter.Domain.Quadrature.LeftEdgey);
p5ValRight = p5(parameter.Domain.Quadrature.RightEdgex,parameter.Domain.Quadrature.RightEdgey);
p5ValUpper = p5(parameter.Domain.Quadrature.UpperEdgex,parameter.Domain.Quadrature.UpperEdgey);
p5ValLower = p5(parameter.Domain.Quadrature.LowerEdgex,parameter.Domain.Quadrature.LowerEdgey);
p5Val = p5(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p5Gradx = 4/dx*p2Val;
p5Grady = 0*p2Val;

p6ValLeft = p6(parameter.Domain.Quadrature.LeftEdgex,parameter.Domain.Quadrature.LeftEdgey);
p6ValRight = p6(parameter.Domain.Quadrature.RightEdgex,parameter.Domain.Quadrature.RightEdgey);
p6ValUpper = p6(parameter.Domain.Quadrature.UpperEdgex,parameter.Domain.Quadrature.UpperEdgey);
p6ValLower = p6(parameter.Domain.Quadrature.LowerEdgex,parameter.Domain.Quadrature.LowerEdgey);
p6Val = p6(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p6Gradx = 0*p6Val;
p6Grady = 4/dy*p3Val;

PVal = [ones(1,length(p6Val));p2Val';p3Val';p4Val';p5Val';p6Val'];
PValLeft = [ones(1,length(p6ValLeft));p2ValLeft;p3ValLeft;p4ValLeft;p5ValLeft;p6ValLeft];
PValRight = [ones(1,length(p6ValRight));p2ValRight;p3ValRight;p4ValRight;p5ValRight;p6ValRight];
PValUpper = [ones(1,length(p6ValUpper));p2ValUpper;p3ValUpper;p4ValUpper;p5ValUpper;p6ValUpper];
PValLower = [ones(1,length(p6ValLower));p2ValLower;p3ValLower;p4ValLower;p5ValLower;p6ValLower];
PGradx = [zeros(1,length(p6Val));p2Gradx';p3Gradx';p4Gradx';p5Gradx';p6Gradx'];
PGrady = [zeros(1,length(p6Val));p2Grady';p3Grady';p4Grady';p5Grady';p6Grady'];


%psi = 1e-4+10*exp(-x.^2*2-y.^2*2)
%psi = 1/4/pi*(1e-10+exp(-10*((x-x0).^2/(2*sigma_x^2)+(y-y0).^2/(2*sigma_y^2))))
%x0 = 0.5;
%y0 = 0.5;
x0 = 0;
y0 = 0;
sigma_x = 0.02;
sigma_y = 0.02;

psi = @(x,y) 1/4/pi*1e-4;

switch parameter.Method.Model
    case 'FullMoment'
        NumberMoments = 1+1/2*parameter.Method.HighestMoment.^2+3/2*parameter.Method.HighestMoment;
        Y = zeros(size(x,1),parameter.Solving.NumberPolynomials,NumberMoments);
        %          Y = cell(NumberMoments,1);
        %          [Y{:}] = deal(zeros(size(x,1),parameter.Solving.NumberPolynomials));
        for i=0:parameter.Method.HighestMoment
            for j=0:i
                k1 = i-j;
                k2 = j;
                ndx = 2+1/2*(i-1).^2+3/2*(i-1)+j;
                M = zeros(size(x,1),parameter.Solving.NumberPolynomials);
                psiMom = @(x,y) 1/2*beta(1/2,1+1/2*(k1+k2))*beta((k1+1)/2,(k2+1)/2)*(1+(-1)^k1+(-1)^(k1+k2)+(-1)^k2)*psi(x,y);
                for k=1:6
                    M(:,k) = (bsxfun(@times,psiMom(QuadraturePointsx,QuadraturePointsy),PVal(k,:)))*Weights/dx/dy;
                end
                %M(:,1) = 1/2*beta(1/2,1+1/2*(k1+k2))*beta((k1+1)/2,(k2+1)/2)*(1+(-1)^k1+(-1)^(k1+k2)+(-1)^k2)/4/pi*max(1e-4,exp(-10*((x(:,1)-x0).^2/(2*sigma_x^2)+(x(:,2)-y0).^2/(2*sigma_y^2)))); %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
                Y(:,:,ndx) = M;
                %                 keyboard
            end
        end
        
    case 'MQM'
        n = parameter.Method.HighestMoment;
        NumberMoments = 1+4.*n+4.*(1/2.*n.^2-1/2.*n);
        i = 1:n;
        A = zeros(1,NumberMoments);
        A(1,1:4*n+1) = [4.*pi,2.*pi./(i+1),(-1).^i.*2.*pi./(i+1),2.*pi./(i+1),(-1).^i.*2.*pi./(i+1)]; %E
        cnt = 0;
        for i=2:n
            for k=1:i-1
                %E
                val = 1/2.*B(1/2,1+i/2).*B((k+1)/2,(i-k+1)/2);
                A(1,4.*n+2+cnt) = val;
                A(1,4.*n+2+nq+cnt) = val.*(-1).^k;
                A(1,4.*n+2+2.*nq+cnt) = val.*(-1).^i;
                A(1,4.*n+2+3.*nq+cnt) = val.*(-1).^(i-k);
                cnt = cnt+1;
            end
        end
        Y = zeros(size(x,1),parameter.Solving.NumberPolynomials,NumberMoments);
        %          Y = cell(NumberMoments,1);
        %          [Y{:}] = deal(zeros(size(x,1),parameter.Solving.NumberPolynomials));
        for i=1:NumberMoments
            Y(:,1,i) = A(i)/4/pi*max(1e-4,exp(-10*((x(:,1)-x0).^2/(2*sigma_x^2)+(x(:,2)-y0).^2/(2*sigma_y^2)))); %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
        end
        
        
end
%keyboard
end

function y = B(u,v)

u(u==0) = 1e-7;
v(v==0) = 1e-7;
y = beta(u,v);
end
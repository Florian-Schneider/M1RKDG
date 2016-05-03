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

p2Val = p2(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p3Val = p3(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p4Val = p4(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p5Val = p5(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);
p6Val = p6(parameter.Domain.Quadrature.ReferencePointsx,parameter.Domain.Quadrature.ReferencePointsy);

PVal = [ones(1,length(p6Val));p2Val';p3Val';p4Val';p5Val';p6Val'];

psi = @(x,y) 1/4/pi*1e-10;


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
                Y(:,:,ndx) = M;
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
            Y(:,1,i) = A(i)*psi; %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
        end
        
        
end
end

function y = B(u,v)

u(u==0) = 1e-7;
v(v==0) = 1e-7;
y = beta(u,v);
end
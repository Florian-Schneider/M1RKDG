function varargout=DirichletBC(t,parameter,~)
persistent argout
%keyboard
if isempty(argout)
    
    x = parameter.Domain.cell_centers(parameter.Domain.n_int+1:end,:);
    Spot2 = 0.45<=x(:,1) & x(:,1) <= 0.55 & x(:,2)<=parameter.Domain.y_lim(1)+parameter.Domain.dy; %Lower Beam
    Spot1 = 0.45<=x(:,2) & x(:,2) <= 0.55 & x(:,1)<=parameter.Domain.x_lim(1)+parameter.Domain.dx; %Left Beam
    
    Psi0 = 100;
    
    
    psi0 = 1/4/pi*1e-10;
    sigma = 0.1/2;
    
    psi1 = @(mue,phi) Psi0.*(exp(-(mue.^2/2/sigma)-(phi.^2/2/sigma))+exp(-(mue.^2/2/sigma)-((phi-2.*pi).^2/2/sigma)));
    psi2 = @(mue,phi) Psi0.*(exp(-(mue.^2/2/sigma)-((phi-pi/2).^2/2/sigma))+exp(-(mue.^2/2/sigma)-((phi-(2.*pi+pi/2)).^2/2/sigma)));
    Omegax = @(mue,phi) sqrt(1-mue.^2).*cos(phi);
    Omegay = @(mue,phi) sqrt(1-mue.^2).*sin(phi);
    switch parameter.Method.Model
        case 'FullMoment'
            NumberMoments = 1+1/2*parameter.Method.HighestMoment.^2+3/2*parameter.Method.HighestMoment;
            argout = cell(NumberMoments,1);
            [argout{:}] = deal(zeros(size(x,1),parameter.Solving.NumberPolynomials));
            for i=0:parameter.Method.HighestMoment
                for j=0:i
                    k1 = i-j;
                    k2 = j;
                    ndx = 2+1/2*(i-1).^2+3/2*(i-1)+j;
                    
                    M = zeros(size(x,1),parameter.Solving.NumberPolynomials);
                    M(:,1) = 1/2*beta(1/2,1+1/2*(k1+k2))*beta((k1+1)/2,(k2+1)/2)*(1+(-1)^k1+(-1)^(k1+k2)+(-1)^k2)*psi0; %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
                    M(Spot1,1) = quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi1(mue,phi),-1,1,0,2*pi,'AbsTol',1e-14);
                    M(Spot2,1) = quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi2(mue,phi),-1,1,0,2*pi,'AbsTol',1e-14);
                    
                    argout{ndx} = M;
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
            
            n = n+1;
            nq = 1/2.*n.^2-1/2.*n;
            i = n;
            A2 = zeros(1,1+4.*n+4.*(1/2.*n.^2-1/2.*n)-NumberMoments);
            A2(1,1:4) = [2.*pi./(i+1),(-1).^i.*2.*pi./(i+1),2.*pi./(i+1),(-1).^i.*2.*pi./(i+1)];
            cnt = 0;
            for i=n:n
                for k=1:i-1
                    %E
                    val = 1/2.*B(1/2,1+i/2).*B((k+1)/2,(i-k+1)/2);
                    A2(1,5+cnt) = val;
                    A2(1,5+nq+cnt) = val.*(-1).^k;
                    A2(1,5+2.*nq+cnt) = val.*(-1).^i;
                    A2(1,5+3.*nq+cnt) = val.*(-1).^(i-k);
                    cnt = cnt+1;
                end
            end
            
            A = [A,A2];
            NumberMoments = 1+4.*n+4.*(1/2.*n.^2-1/2.*n);
            Y = cell(NumberMoments,1);
            [Y{:}] = deal(zeros(size(x,1),parameter.Solving.NumberPolynomials));
            for i=1:NumberMoments
                Y{i}(:,1) = A(i)/4/pi*1e-4; %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
            end
            k1 = 0; k2 = 0; %E
            Y{1}(Spot1,1) = 1;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi1(mue,phi),-1,1,0,2*pi,'AbsTol',1e-14);
            Y{1}(Spot2,1) = 1;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi2(mue,phi),-1,1,0,2*pi,'AbsTol',1e-14);
            k1 = 1; k2 = 0; %Fxp
            Y{2}(Spot1,1) = 1;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi1(mue,phi),-1,1,-pi/2,pi/2,'AbsTol',1e-14);
            Y{2}(Spot2,1) = 0;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi2(mue,phi),-1,1,-pi/2,pi/2,'AbsTol',1e-14);
            k1 = 1; k2 = 0; %Fxm
            Y{3}(Spot1,1) = 0;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi1(mue,phi),-1,1,pi/2,3*pi/2,'AbsTol',1e-14);
            Y{3}(Spot2,1) = 0;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi2(mue,phi),-1,1,pi/2,3*pi/2,'AbsTol',1e-14);
            k1 = 0; k2 = 1; %Fyp
            Y{4}(Spot1,1) = 0;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi1(mue,phi),-1,1,0,pi,'AbsTol',1e-14);
            Y{4}(Spot2,1) = 1;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi2(mue,phi),-1,1,0,pi,'AbsTol',1e-14);
            k1 = 0; k2 = 1; %Fym
            Y{5}(Spot1,1) = 0;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi1(mue,phi),-1,1,pi,2*pi,'AbsTol',1e-14);
            Y{5}(Spot2,1) = 0;%quad2d(@(mue,phi) Omegax(mue,phi).^k1.*Omegay(mue,phi).^k2.*psi2(mue,phi),-1,1,pi,2*pi,'AbsTol',1e-14);
            
            argout = Y;
        case 'SphericalHarmonics'
            NumberMoments = parameter.OptimPara.NN;
            Y = zeros(size(x,1),parameter.Solving.NumberPolynomials,NumberMoments);
            OptimPara = parameter.OptimPara;
            for i=1:NumberMoments
                psiMom = parameter.OptimPara.phiIso(i)*psi0;
                Y(:,1,i) = psiMom;
            end
            psi = psi1(OptimPara.mu,OptimPara.phi)';
            U1 = projectBasis(psi,OptimPara.w',OptimPara.p,3);
            psi = psi2(OptimPara.mu,OptimPara.phi)';
            U2 = projectBasis(psi,OptimPara.w',OptimPara.p,3);
            for i=1:NumberMoments
                Y(Spot1,1,i) = U1(i);
                Y(Spot2,1,i) = U2(i);
            end
            argout = {Y};
    end
end

varargout = argout;

end

function y = B(u,v)

u(u==0) = 1e-7;
v(v==0) = 1e-7;
y = beta(u,v);
end
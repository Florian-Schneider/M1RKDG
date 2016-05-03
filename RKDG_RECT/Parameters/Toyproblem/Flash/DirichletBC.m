function varargout=DirichletBC(t,parameter,~)

x = parameter.Domain.cell_centers(parameter.Domain.n_int+1:end,:);

psi = 1/4/pi*1e-10;

switch parameter.Method.Model
    case 'FullMoment'
        NumberMoments = 1+1/2*parameter.Method.HighestMoment.^2+3/2*parameter.Method.HighestMoment;
        varargout = cell(NumberMoments,1);
        [varargout{:}] = deal(zeros(size(x,1),parameter.Solving.NumberPolynomials));
        for i=0:parameter.Method.HighestMoment
            for j=0:i
                k1 = i-j;
                k2 = j;
                ndx = 2+1/2*(i-1).^2+3/2*(i-1)+j;
                M = zeros(size(x,1),parameter.Solving.NumberPolynomials);
                %                 for m=1:parameter.Solving.NumberPolynomials
                %                     M(:,m) = (m==1)*1/2*beta(1/2,1+1/2*(k1+k2))*beta((k1+1)/2,(k2+1)/2)*(1+(-1)^k1+(-1)^(k1+k2)+(-1)^k2)*(1e-4+19*exp(-x(:,1).^2*2-x(:,2).^2*2)); %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
                %                 end
                M(:,1) = 1/2*beta(1/2,1+1/2*(k1+k2))*beta((k1+1)/2,(k2+1)/2)*(1+(-1)^k1+(-1)^(k1+k2)+(-1)^k2)*psi; %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
                varargout{ndx} = M;
            end
        end
        
    case 'MQM'
        n = parameter.Method.HighestMoment;
        NumberMoments = 1+4.*n+4.*(1/2.*n.^2-1/2.*n);
        nq = 1/2.*n.^2-1/2.*n;
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
            Y{i}(:,1) = A(i)*psi; %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
        end
        varargout = {Y};
        
end

end

function y = B(u,v)

u(u==0) = 1e-7;
v(v==0) = 1e-7;
y = beta(u,v);
end
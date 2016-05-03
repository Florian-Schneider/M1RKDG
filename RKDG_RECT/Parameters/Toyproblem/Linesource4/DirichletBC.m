function varargout=DirichletBC(t,parameter)

x = parameter.Domain.cell_centers(parameter.Domain.n_int+1:end,:);

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
                M(:,1) = 1/2*beta(1/2,1+1/2*(k1+k2))*beta((k1+1)/2,(k2+1)/2)*(1+(-1)^k1+(-1)^(k1+k2)+(-1)^k2)*(1e-4+19*exp(-x(:,1).^2*2-x(:,2).^2*2)); %int(Omegax^k1*Omegay^k2), only cellaverages are prescribed
                varargout{ndx} = M;
            end
        end
        
        
        
end

end
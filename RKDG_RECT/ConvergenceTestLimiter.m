function [] = ConvergenceTestLimiter(k,gamma,nref)
%CONVERGENCETESTLIMITER Summary of this function goes here
%   Detailed explanation goes here

if nargin <3
    nref = 3;
end

nx0 = 5;
ny0 = 5;

U0 = (1-gamma)*[1;1;0]+gamma*[1;0;0]; %(1-gamma)*<b*dirac(Omega-(1;0))> + gamma*(1;0;0)
U1 = 1e-6*((1-gamma)*[1;0;1]+gamma*[1;0;0]); %(1-gamma)*<b*dirac(Omega-(0;1))> + gamma*(1;0;0)

U0 = permute(U0,[2 3 1]);%Vectorized evaluation
U1 = permute(U1,[2 3 1]);
lambda = @(x,y) (cos(2*(x+y)*pi)+1)/2;

U = @(x,y) bsxfun(@times,(1-lambda(x,y)),U0)+bsxfun(@times,lambda(x,y),U1);

L1Error = zeros(length(k),nref+1);
LinfError = L1Error;
thetamax = L1Error;
H =  zeros(1,nref+1);

for c = 1:length(k)
    [PValHQ,PVal,QuadratureHQ,M] = LocalBasis(k(c)); %Local representation of the DG polynomial of degree k+1 on the reference cell
    for i=1:nref+1
        nx = nx0*2^(i-1); %Refine
        ny = ny0*2^(i-1);
        dx = 1/nx;
        dy = 1/ny;
        H(i) = dx;
        [cell_centers_x,cell_centers_y] = meshgrid(dx/2:dx:1-dx/2,dy/2:dy:1-dy/2);
        U_coeff = zeros(nx*ny,3,size(PValHQ,1));
        X = bsxfun(@plus,cell_centers_x(:),(QuadratureHQ.X(:)')*dx);
        Y = bsxfun(@plus,cell_centers_y(:),(QuadratureHQ.Y(:)')*dy);
        Ueval = U(X,Y);
        
        for j=1:size(PValHQ,1) %Project to local basis
            U_coeff(:,:,j) = M(j)^-1*squeeze(sum(bsxfun(@times,Ueval,PValHQ(j,:).*QuadratureHQ.W(:)'),2));
        end
        U_poly_eval = permute(sum(bsxfun(@times,U_coeff,permute(PVal,[3 4 1 2])),3),[1 4 2 3]);
        [U_coeff,theta] = RealizabilityLimiter(U_coeff,U_poly_eval);
        U_poly_eval2 = permute(sum(bsxfun(@times,U_coeff,permute(PValHQ,[3 4 1 2])),3),[1 4 2 3]); %Evaluate at the HQ quadrature for errors
        
        tmp = sqrt(U_poly_eval2(:,:,2).^2+U_poly_eval2(:,:,3).^2)./U_poly_eval2(:,:,1);
        
        if max(tmp(:))>1+10*eps || min(tmp(:)) < 0
            fprintf('%d point-values unrealizable\n',nnz(tmp>1 | tmp < 0)) %Check correctness of limiter
        end
        
        L1Error(c,i) = sum(abs(U_poly_eval2(:,:,1)-Ueval(:,:,1))*QuadratureHQ.W)*dx*dy;
        LinfError(c,i) = max(max(abs(U_poly_eval2(:,:,1)-Ueval(:,:,1))));
        thetamax(c,i) = max(theta);
    end
end

L1Error
LinfError

LinfOrder = bsxfun(@times,-log(LinfError(:,2:end)./LinfError(:,1:end-1)),1./log(H(1:end-1)./H(2:end)))
L1Order = bsxfun(@times,-log(L1Error(:,2:end)./L1Error(:,1:end-1)),1./log(H(1:end-1)./H(2:end)))
thetamax

for c=1:length(k)
data(c).k = k(c);
data(c).Nx = 1./H;
data(c).L1error = L1Error(c,:);
data(c).Linferror = LinfError(c,:);
data(c).nu1 = [NaN,L1Order(c,:)];
data(c).nuinf = [NaN,LinfOrder(c,:)];
data(c).theta = thetamax(c,:);
end
printGrahamsTables(data,'ConvergenceLimiter1.txt')

end

function [PValHQ,PVal,QuadratureHQ,M] = LocalBasis(k)
N = 20; %Number of quadrature nodes per dimension, high order
M = 3; %
p{1} = @(x,y) 1;

if k>1
    p{2} = @(x,y) x;
    p{3} = @(x,y) y;
end
if k>2
    p{4} = @(x,y) p{2}(x,y).*p{3}(x,y);
    p{5} = @(x,y) p{2}(x,y).^2-1/3;
    p{6} = @(x,y) p{3}(x,y).^2-1/3;
end
PValHQ = zeros(length(p),N^2);
PVal = zeros(length(p),M^2);

[z,w]=lgwt(N,-0.5,0.5);
[ZxHQ,ZyHQ] = meshgrid(z);
WHQ = w*w';
WHQ = WHQ(:);
ZxHQ = ZxHQ(:);
ZyHQ = ZyHQ(:);
QuadratureHQ.W = WHQ;
QuadratureHQ.X = ZxHQ;
QuadratureHQ.Y = ZyHQ;

[z,~]=lglnodes(M-1); %Gives nodes on [-1,1]

[Zx,Zy] = meshgrid(z);
Zx = Zx(:);
Zy = Zy(:);


for i=1:length(p)
    PValHQ(i,:) = p{i}(ZxHQ*2,ZyHQ*2); %Polynomials are defined on [-1,1]^2
    PVal(i,:) = p{i}(Zx,Zy); %Polynomials are defined on [-1,1]^2
end
M = sum(bsxfun(@times,WHQ',PValHQ.^2),2);
end

function [U,theta] = RealizabilityLimiter(U,Uval)

U_bar = U(:,:,1); %cell mean

E1 = U_bar(:,1);
Fx1 = U_bar(:,2);
Fy1 = U_bar(:,3);
E2 = Uval(:,:,1);
Fx2 = Uval(:,:,2);
Fy2 = Uval(:,:,3);


tol = 0;

a = (-bsxfun(@plus,E1,-E2).^2*(1-tol) + bsxfun(@plus,Fx1,-Fx2).^2 + bsxfun(@plus,Fy1,-Fy2).^2);
b = 2*(E2.^2*(1-tol) - bsxfun(@times,E1,E2)*(1-tol) - Fx2.^2 + bsxfun(@times,Fx1,Fx2) - Fy2.^2 + bsxfun(@times,Fy1,Fy2));
c = - E2.^2*(1-tol) + Fx2.^2 + Fy2.^2;
q = -1/2*(b+sign(b).*sqrt(b.^2-4*a.*c));
q = real(q);
t1 = q./a;
t2 = c./q;
t1(isnan(t1) | t1>1 | t1<0) = 0;
t2(isnan(t2) | t2>1 | t2<0) = 0;
t1 = max(t1,[],2);
t2 = max(t2,[],2);
theta = max(t1,t2);
%theta = 1;

U(:,:,2:end) = bsxfun(@times,U(:,:,2:end),1-theta);

end

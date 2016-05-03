function [z] = M1(t,Y,parameter)
persistent PValLeft PValRight PValUpper PValLower PGradx PGrady
global PVal
E = Y(:,:,1);
Fx = Y(:,:,2);
Fy = Y(:,:,3);

Eloc = E;
Fxloc = Fx;
Fyloc = Fy;

[EB,FxB,FyB] = feval(parameter.Conditions.BC,t,parameter,Y);

E = [E;EB];
Fx = [Fx;FxB];
Fy = [Fy;FyB];

dx = parameter.Domain.dx;
dy = parameter.Domain.dy;

QuadraturePointsx = parameter.Domain.Quadrature.PhysicalPointsx;
QuadraturePointsy = parameter.Domain.Quadrature.PhysicalPointsy;

EdgeWeightsx = dy/2*[1/6;5/6;5/6;1/6]; %Wrong use lateron; EdgeWeightsx for left/right edges
EdgeWeightsy = dx/2*[1/6;5/6;5/6;1/6];
Weights = parameter.Domain.Quadrature.Weights;


sigma_a = feval(parameter.Physics.sigma_a,QuadraturePointsx,QuadraturePointsy,t);
sigma_s = feval(parameter.Physics.sigma_s,QuadraturePointsx,QuadraturePointsy,t);
Q = feval(parameter.Physics.Source,QuadraturePointsx,QuadraturePointsy,t);

if isempty(PVal) || t==0
    %Basis functions on [-1,1]^2
    %p1 = @(x,y) 1;

  
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

end

LeftNeighbor = parameter.Domain.cell_neighbor_index(:,1);
RightNeighbor = parameter.Domain.cell_neighbor_index(:,2);
LowerNeighbor = parameter.Domain.cell_neighbor_index(:,3);
UpperNeighbor = parameter.Domain.cell_neighbor_index(:,4);

E_h = Eloc*PVal;
E_hLeft = Eloc*PValLeft;
E_hRight = Eloc*PValRight;
E_hLower = Eloc*PValLower;
E_hUpper = Eloc*PValUpper;
E_hLeftN = E(LeftNeighbor,:)*PValRight;
E_hRightN = E(RightNeighbor,:)*PValLeft;
E_hUpperN = E(UpperNeighbor,:)*PValLower;
E_hLowerN = E(LowerNeighbor,:)*PValUpper;

Fx_h = Fxloc*PVal;
Fx_hLeft = Fxloc*PValLeft;
Fx_hRight = Fxloc*PValRight;
Fx_hLower = Fxloc*PValLower;
Fx_hUpper = Fxloc*PValUpper;
Fx_hLeftN = Fx(LeftNeighbor,:)*PValRight;
Fx_hRightN = Fx(RightNeighbor,:)*PValLeft;
Fx_hUpperN = Fx(UpperNeighbor,:)*PValLower;
Fx_hLowerN = Fx(LowerNeighbor,:)*PValUpper;

Fy_h = Fyloc*PVal;
Fy_hLeft = Fyloc*PValLeft;
Fy_hRight = Fyloc*PValRight;
Fy_hLower = Fyloc*PValLower;
Fy_hUpper = Fyloc*PValUpper;
Fy_hLeftN = Fy(LeftNeighbor,:)*PValRight;
Fy_hRightN = Fy(RightNeighbor,:)*PValLeft;
Fy_hUpperN = Fy(UpperNeighbor,:)*PValLower;
Fy_hLowerN = Fy(LowerNeighbor,:)*PValUpper;

[HELeft,HFxLeft,HFyLeft] = H2(E_hLeft,Fx_hLeft,Fy_hLeft,E_hLeftN,Fx_hLeftN,Fy_hLeftN,[-1 0],parameter.Solving.NumberThreads);
[HERight,HFxRight,HFyRight] = H2(E_hRight,Fx_hRight,Fy_hRight,E_hRightN,Fx_hRightN,Fy_hRightN,[1 0],parameter.Solving.NumberThreads);
[HELower,HFxLower,HFyLower] = H2(E_hLower,Fx_hLower,Fy_hLower,E_hLowerN,Fx_hLowerN,Fy_hLowerN,[0 -1],parameter.Solving.NumberThreads);
[HEUpper,HFxUpper,HFyUpper] = H2(E_hUpper,Fx_hUpper,Fy_hUpper,E_hUpperN,Fx_hUpperN,Fy_hUpperN,[0 1],parameter.Solving.NumberThreads);

[FluxX,FluxY] = Flux(E_h,Fx_h,Fy_h,parameter.Solving.NumberThreads);

S = SourceTerm(E_h,Fx_h,Fy_h,sigma_a,sigma_s,Q);

zrho1 = (-(((dx*dy*1)^(-1)*HELeft*EdgeWeightsx+(dx*dy*1)^(-1)*HERight*EdgeWeightsx)+((dx*dy*1)^(-1)*HEUpper*EdgeWeightsy+(dx*dy*1)^(-1)*HELower*EdgeWeightsy))+bsxfun(@times,(dx*dy*1)^(-1)*S{1},PVal(1,:))*Weights);
zrho2 = (-((bsxfun(@times, (dx*dy*1/3)^(-1)*HELeft, PValLeft(2,:))*EdgeWeightsx+bsxfun(@times, (dx*dy*1/3)^(-1)*HERight, PValRight(2,:))*EdgeWeightsx)+(bsxfun(@times, (dx*dy*1/3)^(-1)*HEUpper, PValUpper(2,:))*EdgeWeightsy+bsxfun(@times, (dx*dy*1/3)^(-1)*HELower, PValLower(2,:))*EdgeWeightsy))...
    +(bsxfun(@times, (dx*dy*1/3)^(-1)*FluxX{1}, PGradx(2,:))+bsxfun(@times, (dx*dy*1/3)^(-1)*FluxY{1}, PGrady(2,:)) )*Weights+bsxfun(@times,(dx*dy*1/3)^(-1)*S{1},PVal(2,:))*Weights);
zrho3 = (-((bsxfun(@times, (dx*dy*1/3)^(-1)*HELeft, PValLeft(3,:))*EdgeWeightsx+bsxfun(@times, (dx*dy*1/3)^(-1)*HERight, PValRight(3,:))*EdgeWeightsx)+(bsxfun(@times, (dx*dy*1/3)^(-1)*HEUpper, PValUpper(3,:))*EdgeWeightsy+bsxfun(@times, (dx*dy*1/3)^(-1)*HELower, PValLower(3,:))*EdgeWeightsy))...
    +(bsxfun(@times, (dx*dy*1/3)^(-1)*FluxX{1}, PGradx(3,:))+bsxfun(@times, (dx*dy*1/3)^(-1)*FluxY{1}, PGrady(3,:)) )*Weights+bsxfun(@times,(dx*dy*1/3)^(-1)*S{1},PVal(3,:))*Weights);
zrho4 = (dx*dy*1/9)^(-1)*(-((bsxfun(@times, HELeft, PValLeft(4,:))+bsxfun(@times, HERight, PValRight(4,:)))*EdgeWeightsx+(bsxfun(@times, HEUpper, PValUpper(4,:))+bsxfun(@times, HELower, PValLower(4,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{1}, PGradx(4,:))+bsxfun(@times, FluxY{1}, PGrady(4,:)) )*Weights+bsxfun(@times,S{1},PVal(4,:))*Weights);
zrho5 = (dx*dy*4/45)^(-1)*(-((bsxfun(@times, HELeft, PValLeft(5,:))+bsxfun(@times, HERight, PValRight(5,:)))*EdgeWeightsx+(bsxfun(@times, HEUpper, PValUpper(5,:))+bsxfun(@times, HELower, PValLower(5,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{1}, PGradx(5,:))+bsxfun(@times, FluxY{1}, PGrady(5,:)) )*Weights+bsxfun(@times,S{1},PVal(5,:))*Weights);
zrho6 = (dx*dy*4/45)^(-1)*(-((bsxfun(@times, HELeft, PValLeft(6,:))+bsxfun(@times, HERight, PValRight(6,:)))*EdgeWeightsx+(bsxfun(@times, HEUpper, PValUpper(6,:))+bsxfun(@times, HELower, PValLower(6,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{1}, PGradx(6,:))+bsxfun(@times, FluxY{1}, PGrady(6,:)) )*Weights+bsxfun(@times,S{1},PVal(6,:))*Weights);

zqx1 = (-(((dx*dy*1)^(-1)*HFxLeft*EdgeWeightsx+(dx*dy*1)^(-1)*HFxRight*EdgeWeightsx)+((dx*dy*1)^(-1)*HFxUpper*EdgeWeightsy+(dx*dy*1)^(-1)*HFxLower*EdgeWeightsy))+bsxfun(@times,(dx*dy*1)^(-1)*S{2},PVal(1,:))*Weights);
zqx2 = (-((bsxfun(@times, (dx*dy*1/3)^(-1)*HFxLeft, PValLeft(2,:))*EdgeWeightsx+bsxfun(@times, (dx*dy*1/3)^(-1)*HFxRight, PValRight(2,:))*EdgeWeightsx)+(bsxfun(@times, (dx*dy*1/3)^(-1)*HFxUpper, PValUpper(2,:))*EdgeWeightsy+bsxfun(@times, (dx*dy*1/3)^(-1)*HFxLower, PValLower(2,:))*EdgeWeightsy))...
    +(bsxfun(@times, (dx*dy*1/3)^(-1)*FluxX{2}, PGradx(2,:))+bsxfun(@times, (dx*dy*1/3)^(-1)*FluxY{2}, PGrady(2,:)) +bsxfun(@times,(dx*dy*1/3)^(-1)*S{2},PVal(2,:)))*Weights);
zqx3 = (-((bsxfun(@times, (dx*dy*1/3)^(-1)*HFxLeft, PValLeft(3,:))*EdgeWeightsx+bsxfun(@times, (dx*dy*1/3)^(-1)*HFxRight, PValRight(3,:))*EdgeWeightsx)+(bsxfun(@times, (dx*dy*1/3)^(-1)*HFxUpper, PValUpper(3,:))*EdgeWeightsy+bsxfun(@times, (dx*dy*1/3)^(-1)*HFxLower, PValLower(3,:))*EdgeWeightsy))...
    +(bsxfun(@times, (dx*dy*1/3)^(-1)*FluxX{2}, PGradx(3,:))+bsxfun(@times, (dx*dy*1/3)^(-1)*FluxY{2}, PGrady(3,:)) +bsxfun(@times,(dx*dy*1/3)^(-1)*S{2},PVal(3,:)))*Weights);
zqx4 = (dx*dy*1/9)^(-1)*(-((bsxfun(@times, HFxLeft, PValLeft(4,:))+bsxfun(@times, HFxRight, PValRight(4,:)))*EdgeWeightsx+(bsxfun(@times, HFxUpper, PValUpper(4,:))+bsxfun(@times, HFxLower, PValLower(4,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{2}, PGradx(4,:))+bsxfun(@times, FluxY{2}, PGrady(4,:)) +bsxfun(@times,S{2},PVal(4,:)))*Weights);
zqx5 = (dx*dy*4/45)^(-1)*(-((bsxfun(@times, HFxLeft, PValLeft(5,:))+bsxfun(@times, HFxRight, PValRight(5,:)))*EdgeWeightsx+(bsxfun(@times, HFxUpper, PValUpper(5,:))+bsxfun(@times, HFxLower, PValLower(5,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{2}, PGradx(5,:))+bsxfun(@times, FluxY{2}, PGrady(5,:)) +bsxfun(@times,S{2},PVal(5,:)))*Weights);
zqx6 = (dx*dy*4/45)^(-1)*(-((bsxfun(@times, HFxLeft, PValLeft(6,:))+bsxfun(@times, HFxRight, PValRight(6,:)))*EdgeWeightsx+(bsxfun(@times, HFxUpper, PValUpper(6,:))+bsxfun(@times, HFxLower, PValLower(6,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{2}, PGradx(6,:))+bsxfun(@times, FluxY{2}, PGrady(6,:)) +bsxfun(@times,S{2},PVal(6,:)))*Weights);

zqy1 = (-(((dx*dy*1)^(-1)*HFyLeft*EdgeWeightsx+(dx*dy*1)^(-1)*HFyRight*EdgeWeightsx)+((dx*dy*1)^(-1)*HFyUpper*EdgeWeightsy+(dx*dy*1)^(-1)*HFyLower*EdgeWeightsy))+bsxfun(@times,(dx*dy*1)^(-1)*S{3},PVal(1,:))*Weights);
zqy2 = (-((bsxfun(@times, (dx*dy*1/3)^(-1)*HFyLeft, PValLeft(2,:))*EdgeWeightsx+bsxfun(@times, (dx*dy*1/3)^(-1)*HFyRight, PValRight(2,:))*EdgeWeightsx)+(bsxfun(@times, (dx*dy*1/3)^(-1)*HFyUpper, PValUpper(2,:))*EdgeWeightsy+bsxfun(@times, (dx*dy*1/3)^(-1)*HFyLower, PValLower(2,:))*EdgeWeightsy))...
    +(bsxfun(@times, (dx*dy*1/3)^(-1)*FluxX{3}, PGradx(2,:))+bsxfun(@times, (dx*dy*1/3)^(-1)*FluxY{3}, PGrady(2,:)) +bsxfun(@times,(dx*dy*1/3)^(-1)*S{3},PVal(2,:)))*Weights);
zqy3 = (-((bsxfun(@times, (dx*dy*1/3)^(-1)*HFyLeft, PValLeft(3,:))*EdgeWeightsx+bsxfun(@times, (dx*dy*1/3)^(-1)*HFyRight, PValRight(3,:))*EdgeWeightsx)+(bsxfun(@times, (dx*dy*1/3)^(-1)*HFyUpper, PValUpper(3,:))*EdgeWeightsy+bsxfun(@times, (dx*dy*1/3)^(-1)*HFyLower, PValLower(3,:))*EdgeWeightsy))...
    +(bsxfun(@times, (dx*dy*1/3)^(-1)*FluxX{3}, PGradx(3,:))+bsxfun(@times, (dx*dy*1/3)^(-1)*FluxY{3}, PGrady(3,:)) +bsxfun(@times,(dx*dy*1/3)^(-1)*S{3},PVal(3,:)))*Weights);
zqy4 = (dx*dy*1/9)^(-1)*(-((bsxfun(@times, HFyLeft, PValLeft(4,:))+bsxfun(@times, HFyRight, PValRight(4,:)))*EdgeWeightsx+(bsxfun(@times, HFyUpper, PValUpper(4,:))+bsxfun(@times, HFyLower, PValLower(4,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{3}, PGradx(4,:))+bsxfun(@times, FluxY{3}, PGrady(4,:)) +bsxfun(@times,S{3},PVal(4,:)))*Weights);
zqy5 = (dx*dy*4/45)^(-1)*(-((bsxfun(@times, HFyLeft, PValLeft(5,:))+bsxfun(@times, HFyRight, PValRight(5,:)))*EdgeWeightsx+(bsxfun(@times, HFyUpper, PValUpper(5,:))+bsxfun(@times, HFyLower, PValLower(5,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{3}, PGradx(5,:))+bsxfun(@times, FluxY{3}, PGrady(5,:)) +bsxfun(@times,S{3},PVal(5,:)))*Weights);
zqy6 = (dx*dy*4/45)^(-1)*(-((bsxfun(@times, HFyLeft, PValLeft(6,:))+bsxfun(@times, HFyRight, PValRight(6,:)))*EdgeWeightsx+(bsxfun(@times, HFyUpper, PValUpper(6,:))+bsxfun(@times, HFyLower, PValLower(6,:)))*EdgeWeightsy)...
    +(bsxfun(@times, FluxX{3}, PGradx(6,:))+bsxfun(@times, FluxY{3}, PGrady(6,:)) +bsxfun(@times,S{3},PVal(6,:)))*Weights);


zE = [zrho1,zrho2,zrho3,zrho4,zrho5,zrho6];
zFx = [zqx1,zqx2,zqx3,zqx4,zqx5,zqx6];
zFy = [zqy1,zqy2,zqy3,zqy4,zqy5,zqy6];


z = cat(3,zE,zFx,zFy);

end

function [yE,yFx,yFy] = H2(Eint,Fxint,Fyint,Eext,Fxext,Fyext,nu,N)
%Lax Friedrichs
persistent amax
if isempty(amax)
    amax = 1;
end

[Pxxint,Pxyint,Pyyint,Fxint2,Fyint2] = Closure(Eint,Fxint,Fyint,N);
[Pxxext,Pxyext,Pyyext,Fxext2,Fyext2] = Closure(Eext,Fxext,Fyext,N);


yE = 1/2*((Fxint2+Fxext2)*nu(1)+(Fyint2+Fyext2)*nu(2) -amax*(Eext-Eint)); 
yFx = 1/2*((Pxxint+Pxxext)*nu(1)+(Pxyint+Pxyext)*nu(2)-amax*(Fxext-Fxint)); 
yFy = 1/2*((Pxyint+Pxyext)*nu(1)+(Pyyint+Pyyext)*nu(2)-amax*(Fyext-Fyint)); 



end


function [Pxx,Pxy,Pyy,Fx,Fy,E] = Closure(E,Fx,Fy,N)
tolHack = 1e-10;

Dxx = @(Chi,Nx,Ny) (1-Chi)/2+(3*Chi-1)/2.*Nx.*Nx;
Dxy = @(Chi,Nx,Ny) (3*Chi-1)/2.*Nx.*Ny;
Dyy = @(Chi,Nx,Ny) (1-Chi)/2+(3*Chi-1)/2.*Ny.*Ny;

Fx(E<tolHack) = 0;
Fy(E<tolHack) = 0;
E = max(E,tolHack);

fx = Fx./E;
fy = Fy./E;
absf = sqrt(fx.^2+fy.^2);

fx0 = fx;
fy0 = fy;

fx0(absf>1-tolHack) = fx(absf>1-tolHack)./absf(absf>1-tolHack)*(1-tolHack);
fy0(absf>1-tolHack) = fy(absf>1-tolHack)./absf(absf>1-tolHack)*(1-tolHack);
absf0 = sqrt(fx0.^2+fy0.^2);
Nx = fx0./absf0;
Ny = fy0./absf0;

Nx(absf0<tolHack) = 0;
Ny(absf0<tolHack) = 0;

Nx(isnan(Nx)) = 0;
Ny(isnan(Ny)) = 0;
Nx(isinf(Nx)) = 0;
Ny(isinf(Ny)) = 0;

f = absf0;
chi = 5/3 - (2*(4 - 3*f.^2).^(1/2))/3;

DxyV = Dxy(chi,Nx,Ny).*E;


Pxx = Dxx(chi,Nx,Ny).*E;
Pxy = DxyV;
Pyy = Dyy(chi,Nx,Ny).*E;

end

function [F,G] = Flux(E_h,Fx_h,Fy_h,N)
[Pxx_h,Pxy_h,Pyy_h,Fx_h,Fy_h] = Closure(E_h,Fx_h,Fy_h,N);
F = {Fx_h,Pxx_h,Pxy_h};
G = {Fy_h,Pxy_h,Pyy_h};

end

function [S] = SourceTerm(E_h,Fx_h,Fy_h,sigma_a,sigma_s,Q)
if iscell(Q)
    S = {-sigma_a.*E_h+Q{1},-sigma_a.*Fx_h-2*sigma_s/2.*Fx_h+Q{2},-sigma_a.*Fy_h-2*sigma_s/2.*Fy_h+Q{3}};
else
    S = {-sigma_a.*E_h+4*pi*Q,-sigma_a.*Fx_h-2*sigma_s/2.*Fx_h,-sigma_a.*Fy_h-2*sigma_s/2.*Fy_h};
end
end

function val = p2(x,~)
val = x;
end


function val = p3(~,y)
val = y;
end


function val = p4(x,y)
val = p2(x,y).*p3(x,y);
end


function val = p5(x,y)
val = (p2(x,y).^2-1/3);
end

function val = p6(x,y)
val = (p3(x,y).^2-1/3);
end


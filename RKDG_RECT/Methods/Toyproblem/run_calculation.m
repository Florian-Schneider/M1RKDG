function [fields,val] = run_calculation(parameter)
%RKDG

%Initialization
dt = parameter.Solving.dt;
dt

NumberSteps = ceil(diff(parameter.Solving.tspan)/dt);
NumberPoints = parameter.Domain.n_int;
NumberFrames = 20;%100*ceil(diff(parameter.Solving.tspan)/50);
% NumberFrames = min(NumberFrames,NumberSteps)
%NumberFrames = NumberSteps;

t = linspace(parameter.Solving.tspan(1),parameter.Solving.tspan(2),NumberSteps+1);
dt = t(2)-t(1);
%t = parameter.Solving.tspan(1)+(0:NumberSteps)*dt;
%t(end) = parameter.Solving.tspan(2);
% keyboard
Yn = feval(parameter.Conditions.IC,parameter.Domain.cell_centers(1:NumberPoints,:),parameter);
Yn = Limiter(Yn,parameter,t(1));
Y = cell(1,NumberFrames);
[Y{:}] = deal(Yn);

tframe = 1:NumberFrames;
tframe(1) = t(1);
cnt = 2;
msg = sprintf('t = %2.3f',t(1));
fprintf(msg)
for n=1:NumberSteps+1
    
    if n==NumberSteps+1
        dt = t(end)-t(end-1);
    end
    
    %RK
    Yh0 = Yn;
    %     keyboard;
    [Yh1,~] = Limiter(Yh0 + dt*DifferentialOperator(t(n),Yh0,parameter),parameter,t(n));
    [Yh2,~] = Limiter(3/4*Yh0+1/4*Yh1 + 1/4*dt*DifferentialOperator(t(n),Yh1,parameter),parameter,t(n));
    [Yn,RC] = Limiter(1/3*Yh0+2/3*Yh2 + 2/3*dt*DifferentialOperator(t(n),Yh2,parameter),parameter,t(n));
    
    if mod(n,ceil(NumberSteps/NumberFrames))==0
        Y{cnt} = Yn;
        tframe(cnt) = t(n);
        cnt = cnt+1;
        fprintf(repmat('\b',1,length(msg)));
        msg = sprintf('t = %2.3f',t(n));
        fprintf(msg)
    end
end
Y{cnt} = Yn;
Y = Y(1:cnt);
tframe(cnt) = t(n);
tframe = tframe(1:cnt);
fprintf(repmat('\b',1,length(msg)));
disp('Done');

fields = {'Y','t','tframe','RealizabilityCounter'};
val = {Y,t,tframe,RC};
cd ../..


end

function [Yn] = DifferentialOperator(t,Y,parameter)
Yn = feval(parameter.Method.Space,t,Y,parameter);
end

function [Y,RC]=Limiter(Y,parameter,t)
global PVal
persistent RealCnt
%Flux limiter


if parameter.Method.FluxLimiterFlag~=0
    LeftNeighbor = parameter.Domain.cell_neighbor_index(:,1);
    RightNeighbor = parameter.Domain.cell_neighbor_index(:,2);
    LowerNeighbor = parameter.Domain.cell_neighbor_index(:,3);
    UpperNeighbor = parameter.Domain.cell_neighbor_index(:,4);
    dx = parameter.Domain.dx;
    dy = parameter.Domain.dy;
    E = Y(:,:,1);
    Fx = Y(:,:,2);
    Fy = Y(:,:,3);
    
    [EB,FxB,FyB] = feval(parameter.Conditions.BC,t,parameter,Y);
    
    EB = [E;EB];
    FxB = [Fx;FxB];
    FyB = [Fy;FyB];
    
    ELeftDiff = (E(:,1)-EB(LeftNeighbor,1));
    FxLeftDiff = (Fx(:,1)-FxB(LeftNeighbor,1));
    FyLeftDiff = (Fy(:,1)-FyB(LeftNeighbor,1));
    ERightDiff = -(E(:,1)-EB(RightNeighbor,1));
    FxRightDiff = -(Fx(:,1)-FxB(RightNeighbor,1));
    FyRightDiff = -(Fy(:,1)-FyB(RightNeighbor,1));
    ELowerDiff = (E(:,1)-EB(LowerNeighbor,1));
    FxLowerDiff = (Fx(:,1)-FxB(LowerNeighbor,1));
    FyLowerDiff = (Fy(:,1)-FyB(LowerNeighbor,1));
    EUpperDiff = -(E(:,1)-EB(UpperNeighbor,1));
    FxUpperDiff = -(Fx(:,1)-FxB(UpperNeighbor,1));
    FyUpperDiff = -(Fy(:,1)-FyB(UpperNeighbor,1));
    
    if any(isinf(Y(:)) | isnan(Y(:)))
        keyboard;
    end
    if real(parameter.Method.FluxLimiterFlag)==1
        EMod = E;
        FxMod = Fx;
        FyMod = Fy;
        
        
        
        for j=1:length(E)
            Vx = eye(3);
            Vy = Vx;
            
            ux = [E(j,2);Fx(j,2);Fy(j,2)];
            uL = [ELeftDiff(j,1);FxLeftDiff(j,1);FyLeftDiff(j,1)];
            uR = [ERightDiff(j,1);FxRightDiff(j,1);FyRightDiff(j,1)];
            
            rx = Vx\ux;
            rL = Vx\uL;
            rR = Vx\uR;
            
            rxmod1 = TVBminmod(rx(1),rL(1),rR(1),dx,imag(parameter.Method.FluxLimiterFlag));
            rxmod2 = TVBminmod(rx(2),rL(2),rR(2),dx,imag(parameter.Method.FluxLimiterFlag));
            rxmod3 = TVBminmod(rx(3),rL(3),rR(3),dx,imag(parameter.Method.FluxLimiterFlag));
            
            uy = [E(j,3);Fx(j,3);Fy(j,3)];
            uL = [ELowerDiff(j,1);FxLowerDiff(j,1);FyLowerDiff(j,1)];
            uR = [EUpperDiff(j,1);FxUpperDiff(j,1);FyUpperDiff(j,1)];
            
            ry = Vy\uy;
            rL = Vy\uL;
            rR = Vy\uR;
            
            rymod1 = TVBminmod(ry(1),rL(1),rR(1),dy,imag(parameter.Method.FluxLimiterFlag));
            rymod2 = TVBminmod(ry(2),rL(2),rR(2),dy,imag(parameter.Method.FluxLimiterFlag));
            rymod3 = TVBminmod(ry(3),rL(3),rR(3),dy,imag(parameter.Method.FluxLimiterFlag));
            
            uxmod = Vx*[rxmod1;rxmod2;rxmod3];
            uymod = Vy*[rymod1;rymod2;rymod3];
            EMod(j,2) = uxmod(1);
            EMod(j,3) = uymod(1);
            if abs(EMod(j,2)-E(j,2))>1e-10 || abs(EMod(j,3)-E(j,3))>1e-10
                EMod(j,4:6) = 0;
            end
            FxMod(j,2) = uxmod(2);
            FxMod(j,3) = uymod(2);
            if abs(FxMod(j,2)-Fx(j,2))>1e-10 || abs(FxMod(j,3)-Fx(j,3))>1e-10
                FxMod(j,4:6) = 0;
            end
            
            FyMod(j,2) = uxmod(3);
            FyMod(j,3) = uymod(3);
            if abs(FyMod(j,2)-Fy(j,2))>1e-10 || abs(FyMod(j,3)-Fy(j,3))>1e-10
                FyMod(j,4:6) = 0;
            end
            
        end
        
    else
        EMod = E;
        FxMod = Fx;
        FyMod = Fy;
        
        
        for j=1:length(E)
            tolHack = 1e-5;
            
            e = max(E(j,1),tolHack);
            fx = Fx(j,1)/e;
            fy = Fy(j,1)/e;
            fx(E(j,1)<tolHack) = 0;
            fy(E(j,1)<tolHack) = 0;
            
            n = sqrt(fx^2+fy^2);
            
            HackActive = false;
            if n>1-tolHack
                HackActive = true;
                fx = fx/n*(1-tolHack);
                fy = fy/n*(1-tolHack);
                n = sqrt(fx^2+fy^2);
            end
            if n>tolHack
                nx = fx/n;
                ny = fy/n;
            else
                nx = 0;
                ny = 0;
            end
            f = n;
            
            if f>0
                chi = 5/3 - (2*(4 - 3*f^2)^(1/2))/3;
                dchi = 2*f/sqrt(4-3*f^2);
                a = chi-f*dchi; b = dchi;
                
                dPxxdE = (1-a)/2+(3*a-1)/2*nx^2;
                dPxydE = (3*a-1)/2*nx*ny;
                dPyydE = (1-a)/2+(3*a-1)/2*ny^2;
                dPxxdFx = (-b/2+3*b/2*nx^2+(3*a+3*f*b-1)/f*ny^2)*nx;
                dPxydFx = (-(3*a-1)/2/f*nx^2+(3*a-1+3*f*b)/2/f*ny^2)*ny;
                dPyydFx = (-b/2-(6*a+3*f*b-2)/2/f*ny^2)*nx;
                dPxxdFy = (-b/2-(6*a+3*f*b-2)/2/f*nx^2)*ny;
                dPxydFy = ((3*a-1+3*f*b)/2/f*nx^2-(3*a-1)/2/f*ny^2)*nx;
                dPyydFy = (-b/2+3*b/2*ny^2+(3*a+3*f*b-1)/f*nx^2)*ny;
                
                
                J_x = [0 1 0; dPxxdE dPxxdFx dPxxdFy; dPxydE dPxydFx dPxydFy];
                J_y = [0 0 1; dPxydE dPxydFx dPxydFy; dPyydE dPyydFx dPyydFy];
            else
                J_x = [0 1 0; 1/3 0 0; 0 0 0]; %To avoid numerical instabilities
                J_y = [0 0 1; 0 0 0; 1/3 0 0];
            end
            %
            if ~HackActive
                [Vx,~] = eig(J_x);
                [Vy,~] = eig(J_y);
            else
                Vx = eye(3); %Switch to scalar limiting
                Vy = eye(3);
            end
            
            Vx = real(Vx);
            Vy = real(Vy);
            
            
            ux = [E(j,2);Fx(j,2);Fy(j,2)];
            uL = [ELeftDiff(j,1);FxLeftDiff(j,1);FyLeftDiff(j,1)];
            uR = [ERightDiff(j,1);FxRightDiff(j,1);FyRightDiff(j,1)];
            
            rx = Vx\ux;
            rL = Vx\uL;
            rR = Vx\uR;
            
            rxmod1 = TVBminmod(rx(1),rL(1),rR(1),dx,imag(parameter.Method.FluxLimiterFlag));
            rxmod2 = TVBminmod(rx(2),rL(2),rR(2),dx,imag(parameter.Method.FluxLimiterFlag));
            rxmod3 = TVBminmod(rx(3),rL(3),rR(3),dx,imag(parameter.Method.FluxLimiterFlag));
            
            uy = [E(j,3);Fx(j,3);Fy(j,3)];
            uL = [ELowerDiff(j,1);FxLowerDiff(j,1);FyLowerDiff(j,1)];
            uR = [EUpperDiff(j,1);FxUpperDiff(j,1);FyUpperDiff(j,1)];
            
            ry = Vy\uy;
            rL = Vy\uL;
            rR = Vy\uR;
            
            rymod1 = TVBminmod(ry(1),rL(1),rR(1),dy,imag(parameter.Method.FluxLimiterFlag));
            rymod2 = TVBminmod(ry(2),rL(2),rR(2),dy,imag(parameter.Method.FluxLimiterFlag));
            rymod3 = TVBminmod(ry(3),rL(3),rR(3),dy,imag(parameter.Method.FluxLimiterFlag));
            
            uxmod = Vx*[rxmod1;rxmod2;rxmod3];
            uymod = Vy*[rymod1;rymod2;rymod3];
            EMod(j,2) = uxmod(1);
            EMod(j,3) = uymod(1);
            FxMod(j,2) = uxmod(2);
            FxMod(j,3) = uymod(2);
            
            FyMod(j,2) = uxmod(3);
            FyMod(j,3) = uymod(3);
            
            if norm(ux-uxmod)>1e-10 || norm(uy-uymod)>1e-10
                EMod(j,4:6) = 0;
                FxMod(j,4:6) = 0;
                FyMod(j,4:6) = 0;
            end
        end
    end
    
    Y = cat(3,EMod,FxMod,FyMod);
end
if parameter.Method.RealLimiterFlag~=0
    
    
    if isempty(RealCnt)
        RealCnt = 0;
    end
    
    
    E = Y(:,:,1);
    Fx = Y(:,:,2);
    Fy = Y(:,:,3);
    
    EMean = E(:,1);
    EMean(EMean<1e-14) = 1e-14;
    FxMean = Fx(:,1);
    FyMean = Fy(:,1);
    FMeanNorm = sqrt(FxMean.^2+FyMean.^2);
    FxMean(FMeanNorm>EMean) = FxMean(FMeanNorm>EMean)./FMeanNorm(FMeanNorm>EMean).*EMean(FMeanNorm>EMean)*(1-1e-14);
    FyMean(FMeanNorm>EMean) = FyMean(FMeanNorm>EMean)./FMeanNorm(FMeanNorm>EMean).*EMean(FMeanNorm>EMean)*(1-1e-14);
    E(:,1) = EMean;
    Fx(:,1) = FxMean;
    Fy(:,1) = FyMean;
    E_h = E*PVal;
    Fx_h = Fx*PVal;
    Fy_h = Fy*PVal;
    
    C = sqrt(Fx_h.^2+Fy_h.^2)<=E_h;
    RealCnt = RealCnt + sum(any(C==0,2));
    
    E1 = E(:,1);
    E2 = E_h;
    Fx1 = Fx(:,1);
    Fx2 = Fx_h;
    Fy1 = Fy(:,1);
    Fy2 = Fy_h;
    tol = 0;
    
    a = (-bsxfun(@plus,E1,-E2).^2*(1-tol) + bsxfun(@plus,Fx1,-Fx2).^2 + bsxfun(@plus,Fy1,-Fy2).^2);
    b = 2*(E2.^2*(1-tol) - bsxfun(@times,E1,E2)*(1-tol) - Fx2.^2 + bsxfun(@times,Fx1,Fx2) - Fy2.^2 + bsxfun(@times,Fy1,Fy2));
    c = - E2.^2*(1-tol) + Fx2.^2 + Fy2.^2;
    q = -1/2*(b+sign(b).*sqrt(b.^2-4*a.*c));
    q = real(q);
    t1 = q./a;
    t2 = c./q;
    t1(isnan(t1)) = 1;
    t1(t1>2) = 2;
    t1(t1<-1) = -1;
    t2(isnan(t2)) = 1;
    t2(t2>2) = 2;
    t2(t2<-1) = -1;
    theta = (t1).*(0<=t1 & 1>= t1 & ~(0<=t2 & 1>= t2)) + (t2).*(~(0<=t1 & 1>= t1) & (0<=t2 & 1>= t2)) + (t1).*(0<=t1 & 1>= t1 & (0<=t2 & 1>= t2));
    theta = min(max(theta+tol,[],2),1);
    
    EMod = E;
    FxMod = Fx;
    FyMod = Fy;
    EMod(:,2:6) = bsxfun(@times,(1-theta),E(:,2:6));
    FxMod(:,2:6) = bsxfun(@times,(1-theta),Fx(:,2:6));
    FyMod(:,2:6) = bsxfun(@times,(1-theta),Fy(:,2:6));
    
    
    if parameter.Method.RealLimiterFlag==42
        %            First Order only
        EMod(:,2:end) = 0;
        FxMod(:,2:end) = 0;
        FyMod(:,2:end) = 0;
    end
    
    Y = cat(3,EMod,FxMod,FyMod);
    
    
end

RC = RealCnt;
end

function y = minmod(a1,a2,a3)
index = sign(a1)==sign(a2) & sign(a2)==sign(a3);
y = zeros(size(a1));
y(index) = sign(a1(index)).*min(min(abs(a1(index)),abs(a2(index))),abs(a3(index)));

end

function y = TVBminmod(u,urightdiff,uleftdiff,dx,M)
if nargin < 5
    M = 0;
end
uchange = abs(u)>=M*dx^2;
u(uchange) = minmod(u(uchange),urightdiff(uchange),uleftdiff(uchange));
y = u;
end

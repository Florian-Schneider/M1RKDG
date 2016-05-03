function par = generateParameters(name,GD)
cd(name);

%% Physical Parameters

p.sigma_a = @AbsorptionCoefficient;
p.sigma_s = @ScatteringCoefficient;
p.Source = @Source;

par.Physics = p;

%% Initial Conditions & Boundary Conditions
c.IC = @InitialCondition;
c.BC = @DirichletBC;

par.Conditions = c;

%% Domain

CFL = GD.CFL;
cd('../../..');
d = generateDomain(GD.xspan,GD.yspan,GD.Nx,GD.Ny);
par.Domain = d;

%% Solving
s.epsquer = GD.Time(2)/2;
s.epsmax = 2*s.epsquer;
s.epsS = GD.Time(1);
s.T = 1;
s.tspan = [s.epsS s.epsmax];
%s.dt = CFL*(par.Domain.dx*par.Domain.dy)^2/(par.Domain.dx + par.Domain.dy)^2;
s.dt = CFL*min(par.Domain.dx,par.Domain.dy);
w = 1/12;
dx = par.Domain.dx;
dy = par.Domain.dy;
sigma = max(p.sigma_a(par.Domain.cell_centers(:,1),par.Domain.cell_centers(:,2))+p.sigma_s(par.Domain.cell_centers(:,1),par.Domain.cell_centers(:,2)));
s.dt = min(s.dt,(w*dx*dy)/(dx + dy + sigma*w*dx*dy));
s.NumberPolynomials = 6;

s.TimeStep = s.dt;
s.n = par.Domain.n_int;
s.NumberThreads = 4;


par.Solving = s;
end
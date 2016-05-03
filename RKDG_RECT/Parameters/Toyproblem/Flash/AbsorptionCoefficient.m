function [sigma_a] = AbsorptionCoefficient( x,y, ~ )

R = 1;

n = sqrt(x.^2+y.^2);
sigma_a = 0*(n<=R);
%sigma_a = 0*x+0*y;
end
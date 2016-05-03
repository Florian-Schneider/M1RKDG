function [sigma_a] = AbsorptionCoefficient( x,y, ~ )

O = x>=2 & x<=3 & y>=0 & y<=2;
sigma_a = 50*O;
%sigma_a = 0*x+0*y;
end
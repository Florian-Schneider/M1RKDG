function [sigma_a] = AbsorptionCoefficient( x,y, ~ )
persistent T
if isempty(T)
    %Spot1 = x(:,1)<=3 & x(:,1)>= 2 & x(:,2)<=4 & x(:,2)>= 3;
    %Spot2 = x<=5 & x>= 4 & y<=4 & y>= 3;
    
    Spot1 = x<=2 & x>= 1 & y<=2 & y>=1;
    Spot2 = x<=2 & x>= 1 & y<=4 & y>=3;
    Spot3 = x<=2 & x>= 1 & y<=6 & y>=5;
    Spot4 = x<=6 & x>= 5 & y<=2 & y>=1;
    Spot5 = x<=6 & x>= 5 & y<=4 & y>=3;
    Spot6 = x<=6 & x>= 5 & y<=6 & y>=5;
    Spot7 = x<=3 & x>= 2 & y<=3 & y>=2;
    Spot8 = x<=3 & x>= 2 & y<=5 & y>=4;
    Spot9 = x<=5 & x>= 4 & y<=3 & y>=2;
    Spot10 = x<=5 & x>= 4 & y<=5 & y>=4;
    Spot11 = x<=4 & x>= 3 & y<=2 & y>=1;
    
    
    T = zeros(size(x));
    T(Spot1 | Spot2 | Spot3 | Spot4| Spot5 | Spot6 | Spot7 | Spot8 | Spot9 | Spot10 | Spot11 ) = 10;
    %T(Spot1 | Spot2) = 1;
end

sigma_a = T;
end
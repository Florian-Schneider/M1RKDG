function [Q] = Source( x,y, ~ )

Q = 0*x;
Q(y<=4 & y>= 3 & x<=4 & x>= 3) = 1;
end


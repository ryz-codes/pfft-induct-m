function [M] = analySqRing(R,c)
%ANALYSQRING Analytical self inductance for a square ring with radius R and
%cross sectional thickness c.
%   From Grover, pg 95
%   analySqRing(R,a)
f2 = (c./(2*R)).^2;
P = 4*pi*(0.5*(1+f2/6)*log(8/f2) - 0.84834 + 0.2041*f2);
M = 1e-7.*R.*P;

end


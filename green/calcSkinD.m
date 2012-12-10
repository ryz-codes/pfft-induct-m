function [d] = calcSkinD(f, mu_r, sigma)
%CALCSKIND calculates the skin depth
%   
% Common sigmas:
% Aluminium - 3.5e7
% Copper -  5.8e7
% Iron - 1e7

mu = 1e-7;
d = sqrt(2./(2*pi*f.*mu_r.*mu.*sigma));

end


function [ V, R ] = resist(I,~,L,W,H,IACS)
%RESIST: calculates the PEEC matrix vector product [R]*I for a
%collection of blocks
%   The method simply sets up an internal conductivity and calculates the
%   reistance by volume.
if nargin == 5
    IACS = 1;
end

SIGMA = IACS * 5.8001e7; % International Annealed Copper Standard 100% value

L = sqrt(sum(L.^2,2));
A = sqrt(sum(W.^2,2)) .* sqrt(sum(H.^2,2));

R = L./A/SIGMA;
V = R.*I;


end


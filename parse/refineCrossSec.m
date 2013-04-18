function [Oo,Lo,Wo,Ho] = refineCrossSec(N1,N2,O,L,W,H)
%REFINE_CROSS_SEC Summary of this function goes here
%   Detailed explanation goes here
if nargin == 3 && numel(O) == 4 %if input is a cell
    L = O{2};
    W = O{3};
    H = O{4};
    O = O{1};
end

N = N1*N2;

W = W/N1;
H = H/N2;

Oo = cell(N1,N2);
Lo = cell(N1,N2);
Wo = cell(N1,N2);
Ho = cell(N1,N2);
for ii = 1:N1
    for ij = 1:N2
        Oo{ii,ij} = O + (ii-1)*W + (ij-1)*H;
        Lo{ii,ij} = L;
        Wo{ii,ij} = W;
        Ho{ii,ij} = H;
    end
end

Lo = repmat(L,N,1);
Wo = repmat(W,N,1);
Ho = repmat(H,N,1);
Oo = vertcat(Oo{:});

if nargout == 1
    Oo = {Oo, Lo, Wo, Ho};
end
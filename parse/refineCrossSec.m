function [Oo,Lo,Wo,Ho] = refineCrossSec(N1,N2,O,L,W,H)
%REFINE_CROSS_SEC Refines rectangular filaments into smaller sub filaments. 
% Allows a fill factor to be entered in order for the effects of insulation
% to be simulated.
if nargin == 3 && numel(O) == 4 %if input is a cell
    L = O{2};
    W = O{3};
    H = O{4};
    O = O{1};
end
if nargin == 0 % demo
    N1 = 4; N2 = 4;
    O = [0 0 0];
    L = [10 0 0];
    W = [0 1 0];
    H = [0 0 0.5];
end

% Fill factor is hard-set.
% wdiam = 0.8128e-3;
% (wdiam/2)^2*pi*16/wid/hei
fillfac = 0.6924;

% Total number of new filaments
N = N1*N2;

% Width reduction due to refinement
W = W/N1;
H = H/N2;

% Width reduction due to fill factor
wfac = sqrt(fillfac); % conductor thickness wrt the width of fill filament
dw = 0.5-wfac/2; %insulation thickness wrt to width of full filament

Oo = cell(N1,N2);
for ii = 1:N1
    for ij = 1:N2
        Oo{ii,ij} = O + (ii-1+dw)*W + (ij-1+dw)*H;
    end
end
Oo = vertcat(Oo{:});
Lo = repmat(L,N,1);
Wo = repmat(W*wfac,N,1);
Ho = repmat(H*wfac,N,1);


if nargout == 1
    Oo = {Oo, Lo, Wo, Ho};
end
if nargout == 0
   showFils(Oo,Lo,Wo,Ho); 
end
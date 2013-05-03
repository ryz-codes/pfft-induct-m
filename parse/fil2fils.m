function [O,L,W,H] = fil2fils(Oin, Lin, Win, Hin, pxrad)
%Fil2FILS Like seg2fils, but divides a single filament into subfilaments that can
%be treated with eddy current treatment
persistent pxrad2
persistent gi
persistent gj
persistent gw
persistent gh

if nargin == 0
    Oin = [0,0,0]; Lin = [10,0,0]; Win = [0,2,0]; Hin = [0,0,3];
    pxrad = 4;
end

if pxrad ==1
    O = Oin; L = Lin; W = Win; H = Hin;
    return;
end

% Pixellate the cross section
if ~isequal(pxrad,pxrad2)
    pxrad2 = pxrad;
    %tmp = 0.5*(sin(linspace(-pi/2,pi/2,2*pxrad-1))+1);
    tmp = linspace(-pxrad/2,pxrad/2,2*pxrad-3);
    tmp = [0,tanh(tmp)/2+0.5,1];
    [gi, gj] = ndgrid(tmp(1:end-1));
    [gw, gh] = ndgrid(diff(tmp));
end
N = numel(gi);

O = bsxfun(@plus,Oin,gi(:)*Win + gj(:)*Hin);
L = repmat(Lin, N, 1);
W = gw(:)*Win;
H = gh(:)*Hin;

if nargout ==0 
showFils(O,L,W,H);
size(O,1)
axis off
end

end


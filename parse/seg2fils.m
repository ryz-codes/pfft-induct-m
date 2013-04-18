function [O,L,W,H] = seg2fils(from, to, rad, pxrad)
%SEG2FILS Summary of this function goes here
%   Detailed explanation goes here
%from = [0 0 0]; to = [1 0 0];
%rad = 0.3;
%pxrad = 5;
persistent pxrad2
persistent gi
persistent gj
dr = rad/pxrad;

ux = (to-from); % lengthwise along the path
ux = ux./norm(ux); 
uy = ux*[0 1 0; -1 0 0; 0 0 1]; 
uz = cross(ux,uy);

% Pixellate the cross section
if ~isequal(pxrad^2,pxrad2)
    pxrad2 = pxrad^2;
    [gi, gj] = ndgrid(-pxrad:pxrad-1);
    idx = (gi.^2+gj.^2 >= pxrad2);
    gi(idx) = []; gj(idx) = [];
end

O = gi(:)*uy + gj(:)*uz;
O = bsxfun(@plus,from+dr*(uy/2-uz/2),dr*O);
N = size(O,1);
L = repmat(to-from,N,1);
W = repmat((dr*uy),N,1);
H = repmat((dr*uz),N,1);

%showFils(O,L,W,H);

end


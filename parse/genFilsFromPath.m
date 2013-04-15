function [O,L,W,H] = genFilsFromPath(path, wid, hei)
%GENFILSFROMPATH Summary of this function goes here
%   Detailed explanation goes here
%
% Path should be a 3-column matrix
%
% enter wid and hei as scalars to set them horizontal and vertical
% 
% (not implemented)
% wid and hei should be a 2-vector, in the plane of the cross section of a 
% unit length vector.


N_tot = size(path,1)-1;

O = path(1:N_tot,:);
L = path(2:N_tot+1,:) - O;
len = sqrt(sum(L.^2,2));
if nargin == 1 % width & height not provided.
    wid = mean(len)*1e-3; % set aspect ratio to 1000:1
    hei = mean(len)*1e-3;
end

uL = bsxfun(@rdivide,L, len); % unit vectors
% MAKE WIDTH
% width is always parallel to the xy plane
% 1. project to xy plane
% 2. rotate along the z axis
uLp = bsxfun(@rdivide,uL(:,1:2),sqrt(sum(uL(:,1:2).^2,2))); % proj to xy plane
W = uLp * [0 1; -1 0]; % rotate 
W = [W zeros(size(W,1),1)];


% MAKE HEIGHT
% height is rotation with W as axis.
% 1. cross product W with uL
H = cross(uL,W,2); % rotate with y

% dimensionalize.
W = W *wid;
H = H *hei;

% Shift everything down W and H
O = O - 0.5*W - 0.5*H;

if nargout == 1
    O = {O,L,W,H};
end

end


function [litz_verts] = path2litz(verts, rad, num, len)
%PATH2LITZ Converts a single vertex path into a group of litz paths
%   rad = radius of the litzing, num = number of litz paths, 
%   len = number of vertices per full litz rotation.
%   The litzing pattern is r = cos(theta/4). Where theta comes from inside
%   the widing.
if nargin ==0
verts = genSpiralPath([0.5,1,1.5],0,100);
showPath(verts);
len = 10;
num = 10;
rad = 0.1;
end

Nvert = size(verts,1);

% Make an array of angles
ang = bsxfun(@plus,(1:Nvert).',linspace(0,len*4,num+1));

% Delete last column (repeat of first) and scale. Allocate litzing angle
% and radius.
ang = ang(:,1:end-1)/len*2*pi;
r = rad*cos(ang/4);

% Figure out the directions (uy is along width, uz is along height)
L = verts(2:end,:) - verts(1:end-1,:); L = [L; L(1,:)];
W = L * [0 1 0; -1 0 0; 0 0 1];
H = cross(W,L,2);
uy = bsxfun(@rdivide,W,sqrt(sum(W.^2,2)));
uz = bsxfun(@rdivide,H,sqrt(sum(H.^2,2)));

% Create the displacement vectors for each litz.
litz_verts = cell(1,num);
for ii = 1:num
    litz_verts{ii} = verts + bsxfun(@times,uy,r(:,ii).*cos(ang(:,ii))) ...
        + bsxfun(@times,uz,r(:,ii).*sin(ang(:,ii)));
end

% Display all wires
if nargin == 0
showPath(vertcat(litz_verts{:}));
axis equal
end
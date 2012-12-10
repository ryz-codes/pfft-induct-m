function [O L W H] = genCircFils( R, Z, N, wid, hei )
%genCircFils Generates a circle of filaments
%   Detailed explanation goes here

O = cell(numel(R),1);
L = cell(numel(R),1);
W = cell(numel(R),1);
H = cell(numel(R),1);

for ii = 1:numel(R)
    ang = linspace(0,2*pi,N+1);
    verts = [R(ii)*cos(ang); R(ii)*sin(ang); ones(size(ang))*Z].';

    O{ii} = verts(1:N,:);
    L{ii} = verts(2:N+1,:) - O{ii};
    ul = L{ii} ./ (sqrt(sum(L{ii}.^2,2)*[1 1 1]));
    W{ii} = ul * [0 -1 0; 1 0 0; 0 0 1] * wid; % rotate
    H{ii} = repmat([0 0 1],N,1) * hei; % copy
end

O = vertcat(O{:});
L = vertcat(L{:});
W = vertcat(W{:});
H = vertcat(H{:});
    
% Shift everything down W and H
O = O - 0.5*W - 0.5*H;

end


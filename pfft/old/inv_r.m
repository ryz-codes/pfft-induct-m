%INV_R vectorized 1/r green's function, with zero at the origin.
% Please input the arguments as a cell, one cell for each dimension.
% 
% Has an internal setting on the minimum distance before setting things to
% zero.
function g = inv_r(X)
MIN_DIST = 1e-1;
Z_DIST_SQ = 1e-6;

assert(iscell(X));

dist = Z_DIST_SQ;
for ii = 1:numel(X)
    dist = dist + X{ii}.^2;
end
dist = sqrt(dist);
g = zeros(size(dist));

% Find the non-zero points
z_idx = dist > MIN_DIST;
g(z_idx) = 1./dist(z_idx);
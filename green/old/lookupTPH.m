function [G] = lookupTPH(ri,zh,zt,H,T,r,z,zp)
%LOOKUPTPH Summary of this function goes here
%   Table details
%   r,z,zp - matrix of r, z, zp of the points to be examined.
%   H - Hankel table, spanned by ri along its columns and zh along its
%       rows.
%   T - Toeplitz table, spanned by ri along its columns and zt along its
%       rows.
%   ri - r values along the columns of the table
%   zh - hankel row values
%   zt - toeplitz row values
xh = z+zp;
xt = abs(z-zp);

if numel(xh) == 1
    xh = xh * ones(size(r));
    xt = xt * ones(size(r));
end


Gh = interp2(ri,zh,H,r,xh);
Gt = interp2(ri,zt,T,r,xt);
% Gh = interp2(xh,r,H,zh,ri);
% Gt = interp2(xt,r,T,zt,ri);

G = Gh+Gt;
end


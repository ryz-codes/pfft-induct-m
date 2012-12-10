function [G T H] = makeTPH(r,z,zp)
% MAKETPH makes the Toeplitz plus Hankel green's function matrix with the
% Laplacian Green's function

% Define physical system
% interface is located at z=d
% ==Above the surface
% charge 1 appears at zp
% charge 2 appears at -zp+2*d as (mu-1)/(mu+1)
% ==Below the surface
% charge 1 appears at zp as 2mu/(mu+1)
d = 0;
mu = 30;


assert(all(length(zp(:))==1),'zp must contain one element'); 
N = length(r); %number of columns to fetch
M = length(z); %number of rows

% make matrices
r = kron(r(:).',ones(M,1)); % copy row vector
z = kron(ones(1,N),z(:)); % copy column vector
zp = kron(zp,ones(M,N)); % copy column vector

% above surface
Tu = gf(r,z-zp).*(z>=d);
Hu = gf(r,z+zp-2*d)*(mu-1)/(mu+1).*(z>=d);
Gu = Tu+Hu;

% below surface
Tb = gf(r,z-zp)*2*mu/(mu+1).*(z<d);
Gb = Tb;

% sum
T = Tu+Tb;
H = Hu;
G = Gu+Gb;

    function out = gf(r,z) %point green's function
        out = 1./sqrt(r.^2 + z.^2);
    end


end
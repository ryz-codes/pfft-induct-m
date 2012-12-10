function [A r] = blk( rN,rbnd,zN,dz,psn )
%BLK Summary of this function goes here
%   Detailed explanation goes here
assert(nargin>=4,'Cannot make block without rN, rbnd, zN and zbnd');
if nargin == 4
    psn = 0;
end

% rN = 60;
% rbnd = 1; % r upper bound.
% zN = 600;
% dz = 2.5e-4;
% psn = -1j;

%%
%--------------------------------------------------------------------------
% R BLOCK
%
% Generate one block of rN x rN matrix for d2/dr2 + 1/r d/dr A = 0 
% with non-uniform, polynomially spaced grid.
%
% Output - Ar, the single bock structure.
%--------------------------------------------------------------------------
x = linspace(0,1,rN).'; % Define a uniformly space variable x.
dx = x(2)-x(1); % Extract x spacing.
ord = 3; % Define order of polynomial relationship between x and r
r = (x.^ord + 0.5./(rN-1).^ord) * rbnd; % Shift r to enforce Neumann
dr = r(2)-r(1); % Extract r spacing.

% 1 x rN row matrices of the entries involved.
% r_fwd and r_rev to deal with the second derivative
% r_dif to deal with the first derivative
r_fwd = 1./(ord*(x+dx/2).^(ord-1).*(ord*x.^(ord-1)).*dx^2.*rbnd^2);
r_rev = 1./(ord*(x-dx/2).^(ord-1).*(ord*x.^(ord-1)).*dx^2.*rbnd^2);
r_dif = 1./(2*rbnd*dx .* r .* (ord*x.^(ord-1)));

% Diagonals written in a full flat matrix
B = [r_fwd  +  r_dif, ...
     r_rev  -  r_dif, ...
     -r_fwd -  r_rev];

% Form (rN-1) x rN matrix without x=end+1 row to implicitly enforce Dirichlet.
Ar = spdiags(B(2:end,:), ...
    [2 0 1], rN-1,rN);

% Enfore Neuman condition at x=0. Construct the matrix without
% that row, then ENFORCE: x(0) = x(1).
bnd_factor = 1/dr^2+1/(2*dr*r(1));
neu_bnd = sparse([1; 1],[1; 2],...
    [-1; 1]*bnd_factor,...
    1,rN);
Ar = [neu_bnd;Ar];

%%
%--------------------------------------------------------------------------
% Z BLOCK
%
% Generates d2/dz2 A - jw mu sig A = 0 structure
%
% Output - Az, the single bock structure for z.
%--------------------------------------------------------------------------
z_fwdrev = 1./ dz.^2;
B = ones(zN,3) * z_fwdrev;
B(:,2) = -2 * B(:,2);
Az = spdiags(B,[-1 0 1],zN,zN);

%%
%--------------------------------------------------------------------------
% Full BLOCK
%
% Generates Ar + Az + psn A = 0 structure
%
% Input - psn, the Poisson term.
% Output - A, the total block structure for z.
%--------------------------------------------------------------------------
A = kron(speye(zN),Ar)+kron(Az,speye(rN));
if psn ~= 0
    A = A + psn*speye(rN*zN);
end

end


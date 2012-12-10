function [Aout,r,zout] = edd_fast(z_range,zp)
% This FDM code solves the following system:
%       1/r dA/dr + d^2A/d^2r + d^2A/d^z - jw mu sigma A = -mu J
% 
% Boundary conditions are explictly enforced for both continuity of A and
% continuity of dA/dn. The excitation is achieved through boundary
% conditions rather than domain functions.
%
% This version conditions the matrices to be relatively well conditioned.
%
% Non-uniform FDM techniques were used according to 
%   Eugenia Kálnay de Rivas, "On the use of nonuniform grids in 
%   finite-difference equations", Journal of Computational Physics, 
%   Volume 10, Issue 2, October 1972, Pages 202-210
% with a polynomial grid structure. 
% 
% Adapted from the folder "0623-nonuniform"



if nargin == 1
    z_range = zp;
end

tic
%% Setup
% Define region
mu = 4*pi*1e-7;
w = 2*pi*2e4;
mu_r = 130;
sig = 8e6;
d=abs(zp);
%d = 5e-3;
skin_d = sqrt(2/w/sig/mu_r/mu);

d

% Non-uniform grid polynomial order
ord = 3;

% Easy setup
MULTIPLIER = 300;
ITER = 0;
disp('  ');disp('  ');

% Definition the FDM grid
rN = 60;
r_max = 1;

% multiresolution scheme.
% freespace 
zN1 = 6*MULTIPLIER;
z_min1 = -0.5;
z_max1 = 0; % the coarse grid ends at d.
% material
zN2 = 2*MULTIPLIER;
z_min2 = 0; % the fine grid starts at d.
if ~isinf(skin_d)
    z_max2 = skin_d*5; % goes for 5 skin depths.
else
    z_max2= -z_min1;
end

zN = zN1+zN2;

% Make r and z grids
z1 = linspace(z_min1,z_max1,zN1);
z_step1 = (z_max1-z_min1)/(zN1-1);

z2 = linspace(z_min2,z_max2,zN2);
z_step2 = (z_max2-z_min2)/(zN2-1);

% make uniformly sampled space
x = linspace(0,1,rN);
x_step = x(2)-x(1);
x_all = kron(ones(1,zN),x);

% make non-uniformly sampled space
r = (x.^ord + 0.5./(rN-1).^ord) * r_max;
r_all = kron(ones(1,zN),r);
z_all = [kron(z1,ones(1,rN)), kron(z2,ones(1,rN))];

% Timing
disp(['Setup: ' mat2str(toc)]); tic

%% Fill matrix
% Figure out which index number refers to which cell
% Boundaries
n_top=1:rN; % top
n_bot = (zN-1)*rN + (1:rN); % bottom
n_lef1 = (rN+1):rN:((zN1-2)*rN+1); % left boundary in air
n_lef2 = ((zN1+1)*rN+1):rN:((zN-2)*rN+1); % left in material
n_rig = rN*2:rN:(zN-1)*rN; % right (without top and bot rows)

% The eddy current boundary is the last row of the freespace domain and the
% first row of the material domain. 
n_int = (zN1-1)*rN + (1:rN-1); 
n_int2 = (zN1)*rN + (1:rN-1);

% Freespace
n_free = ones(zN1-2,1) * (2:rN-1) + ((1:zN1-2)*rN).' * ones(1,rN-2);
n_free = n_free(:).';

% Eddy current region
n_eddy = ones(zN2-2,1) * (2:rN-1) + ((zN1+1:zN-2)*rN).' * ones(1,rN-2);
n_eddy = n_eddy(:).';

%%

% Create the free-space domain entries
z_mult = 1./z_step1^2*ones(1,length(n_free));
r_mult_fwd = 1./(ord*(x_all(n_free)+x_step/2).^(ord-1).*(ord*x_all(n_free).^(ord-1)).*x_step^2.*r_max^2);
r_mult_rev = 1./(ord*(x_all(n_free)-x_step/2).^(ord-1).*(ord*x_all(n_free).^(ord-1)).*x_step^2.*r_max^2);
r_dif = 1./(2*r_max*x_step .* r_all(n_free) .* (ord*x_all(n_free).^(ord-1)));

rows1 = kron([1 1 1 1 1],n_free);
cols1 = [n_free+1,n_free-1,n_free+rN,n_free-rN,n_free];
ent1 = [r_mult_fwd+r_dif,r_mult_rev-r_dif,z_mult,z_mult,-r_mult_fwd-r_mult_rev-2*z_mult];

% Create the eddy domain entries
z_mult = ones(1,length(n_eddy));
r_mult_fwd = z_step2^2./(ord*(x_all(n_eddy)+x_step/2).^(ord-1).*(ord*x_all(n_eddy).^(ord-1)).*x_step^2.*r_max^2);
r_mult_rev = z_step2^2./(ord*(x_all(n_eddy)-x_step/2).^(ord-1).*(ord*x_all(n_eddy).^(ord-1)).*x_step^2.*r_max^2);
r_dif = z_step2^2./(2*r_max*x_step .* r_all(n_eddy) .* (ord*x_all(n_eddy).^(ord-1)));

rows2 = kron([1 1 1 1 1],n_eddy);
cols2 = [n_eddy+1,n_eddy-1,n_eddy+rN,n_eddy-rN,n_eddy];
ent2 = [r_mult_fwd+r_dif,r_mult_rev-r_dif,z_mult,z_mult, ...
    -r_mult_rev-r_mult_fwd-2*z_mult-z_step2^2*1j*w*mu*mu_r*sig]; % Helmholtz term

% Create the matrix indices for the outer boundaries (Dirichlet)
rows3 = [n_top, n_bot, n_rig];
cols3 = [n_top, n_bot, n_rig];
ent3 = ones(size(rows3)); 

% calculate r_step at origin
r_step = r(2)-r(1);

% Inner boundary in air (Neumann)
rows6 = [n_lef1 n_lef1 n_lef1 n_lef1];
cols6 = [n_lef1-rN n_lef1+rN n_lef1+1 n_lef1];
r_dif = z_step1^2./(1*r_step * r_all(n_lef1));
z_mult = ones(1,length(n_lef1));
r_mult_fwd = z_step1^2./r_step^2.*ones(1,length(n_lef1));
ent6 = [z_mult,z_mult,r_mult_fwd+r_dif,-r_mult_fwd-2*z_mult-r_dif];

% Inner boundary in material (Neumann)
rows7 = [n_lef2 n_lef2 n_lef2 n_lef2];
cols7 = [n_lef2-rN n_lef2+rN n_lef2+1 n_lef2];
r_dif = z_step2^2./(1*r_step * r_all(n_lef2));
z_mult = ones(1,length(n_lef2));
r_mult_fwd = z_step2^2./r_step^2.*ones(1,length(n_lef2));
ent7 = [z_mult,z_mult,r_mult_fwd+r_dif, ...
    -r_mult_fwd-2*z_mult-r_dif-1j*z_step2^2*w*mu*mu_r*sig];

% Continuity of 1/mu dA/dn at the boundary.
rows4 = [n_int,n_int,n_int,n_int];
cols4 = [n_int,n_int-rN,n_int+2*rN,n_int+rN];
temp1 = ones(1,length(n_int));
temp2 = z_step1./z_step2./mu_r.*ones(1,length(n_int));
ent4 = [temp1,-temp1,-temp2,temp2]; 

% Continuity of A at the boundary
rows5 = [n_int2,n_int2];
cols5 = [n_int2-rN,n_int2];
ent5 = [-ones(1,length(n_int2)), ones(1,length(n_int2))];

% Piece together the sparse matrix format
rows = [rows1,rows2,rows3,rows4,rows5,rows6,rows7];
cols = [cols1,cols2,cols3,cols4,cols5,cols6,cols7];
entries = [ent1,ent2,ent3,ent4,ent5,ent6,ent7];
H = sparse(rows,cols,entries,rN*zN,rN*zN);

% Timing
disp(['Matrix fill: ' mat2str(toc)]); tic

%% Fill RHS
b = zeros(rN*zN,1);
b(n_int) = z_step1*mu./4./pi.*(z_all(n_int)+d)./(sqrt(r_all(n_int).^2 + (z_all(n_int)+d).^2).^3);
b(n_int2) = mu./4./pi./sqrt(r_all(n_int2).^2 + (z_all(n_int2)+d).^2);

% Timing
disp(['Column fill: ' mat2str(toc)]);

%% Clean up
clear rows cols entries rows1 rows2 cols1 cols2 ent1 ent2 n_top n_bot 
clear n_lef n_rig n_free z_mult r_mult r_diff r_all z_all n_eddy


% if ITER ~= 1
    %% Back-sub Solve
    tic
    A = H \ b;
    % Timing
    disp(['A \ b time: ' mat2str(toc)]);
% else
%     %% Iterative solve
%     
%     tol=1e-3;
%     
%     %Precond
%     tic; [L,U] = luinc(H,struct('droptol',0.1,'milu',1,'thresh',1));  
%     disp(['Preconditioner: ' mat2str(toc)]); tic
%     
%     % GMRES
%     tic
%     [A,flag,relres,iter,resvec] = gmres(H,b,[],tol,50,L,U);
%     disp(['GMRES: ' mat2str(iter(2)) ' iterations in ' mat2str(toc)]);
%     figure(3);semilogy(0:iter(2),resvec./resvec(1));
%     if flag==1
%         warning('GMRES: Did not converge.');
%     end
% end

A = reshape(A,rN,zN); % morph the solution back into readable form

% Extract data
z = [z1 z2];
ind1 = find(z>=z_range(1), 1, 'first');
if numel(z_range)==1
    ind2 = ind1;
else
    ind2 = find(z<=z_range(2), 1, 'last');
end
Aout = A(:,ind1:ind2);
zout = z(ind1:ind2);
%% Extract data
% ind1 = find(z1>=-d, 1, 'first');
% A1 = A(:,ind1);
% d1 = z1(ind1);
% A2 = A(:,find(z1==0, 1, 'first'));
% save('FDM_data','A1','A2','A','r','z1','R','d1','d','w');
% figure(1);plot(r,real(A1),r,imag(A1));
% draw_response_fht
% figure(3)


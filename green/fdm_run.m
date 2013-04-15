function [xout rout zout] = fdm_run(z_range,zp,L)
%FDM_RUN - Runs FDM to retrieve the Green's function for one particular zp
% This FDM code solves the following system:
%       1/r dA/dr + d2A/dr2 + d2A/dz2 - jw mu sigma A = -mu J
%
% With the boundary conditions:
%        A2           -  A1          = 0
%        1/mu2 dA2/dz - 1/mu1 dA1/dz = 0
%
% For the input of
%       J = delta(x) delta(y) delta(z-zp) 
%
% This excitation is achieved through boundary conditions rather than 
% domain functions.
%
% Non-uniform FDM techniques were used according to 
%   Eugenia Kálnay de Rivas, "On the use of nonuniform grids in 
%   finite-difference equations", Journal of Computational Physics, 
%   Volume 10, Issue 2, October 1972, Pages 202-210
% with a polynomial grid structure. 
% 
% Inputs:
%    L    - the definitions of the multilayered structure (see defaultL for 
%           descriptions of fields)
%
% Outputs:
%    g    - A green's function table object. The table can be accessed by
%           typing the following command:
%                   G = g.lookup(r,z,zp);
% 
% Other m-files required: fdm_run (the fdm algorithm that generates the table)
% Subfunctions: none
% MAT-files required: none
%
% See also: defaultL.m
%
%
% Author: Richard Y Zhang
% Massachusetts Institute of Technology
% Email: ryz@mit.edu  
% Oct 2012; Last revision: 8th Oct 2012


%% SYSTEM DEFAULTS
% Definition the radial FDM grid
rN = 60;
rbnd = 1;

% Non-uniform grid polynomial order
ord = 3;

%% IMPORT DATA FROM STRUCT
% Read "defaultL.m" for details on these variables
w = L.w;               mu_r = L.mu_r;     sig = L.sig;                   
mu = L.mu;             zN = L.zN;         coil_layer = L.coil_layer;     
layerN = L.layerN;     bnds = L.bnds; 

%% MAIN FUNCTION
% Make a vector of z values and boundaries
z_val = [];
bnd_ind = zeros(layerN,2); % index of the upper and lower boundaries
for ii = 1:layerN
    temp = linspace(bnds(ii),bnds(ii+1),zN(ii));
    dz(ii) = temp(2) - temp(1); % figure out the spacing
    
    bnd_ind(ii,1) = length(z_val)+1; % first index in this layer
    bnd_ind(ii,2) = length(z_val)+zN(ii); % last index in this layer
    z_val = [z_val, temp];
end

% Make A block for each layer region
Alay = cell(layerN,1);
psn = -1j*w*mu.*sig.*mu_r; % poisson factors

[Ar r] = blk_ar; % radial blocks
if layerN >2
    for ii = 2:(layerN-1)
        Alay{ii} = blk_asm(Ar,zN(ii)-2,dz(ii),psn(ii));
        % delete two z elements to make space for one boundary on each side
    end
end
for ii = [1,layerN]
    Alay{ii} = blk_asm(Ar,zN(ii)-1,dz(ii),psn(ii));
    % delete one z element to make space for the boundary.
end

% Connect the boundaries
fac = 1./dz./mu_r; % boundary factors
A = Alay{1};
for ii = 2:layerN
    A = blk_bnd(A,Alay{ii},fac([ii-1,ii]));
end

% DEBUG, make sure that z_val and A match or else all index extraction are
% wrong
assert((length(A)/rN) == length(z_val));

% Find the excitation boundaries
bnd_above = bnds(coil_layer+1);
bnd_below = bnds(coil_layer);
bnd_a_ind = rN*(bnd_ind(coil_layer,2)-1);
bnd_b_ind = rN*(bnd_ind(coil_layer,1)-2);

% Logic for doing boundary above or below
if coil_layer == 1 %|| ... Coils is at the bottom or...
        %(mu_r(coil_layer-1)==1 && sig(coil_layer-1)==0) 
        % next layer is air (assume air forever)
    bnd_below = [];
elseif coil_layer == layerN% || ... Coil is at the top or...
        %(mu_r(coil_layer+1)==1 && sig(coil_layer+1)==0) 
        % next layer is air (assume air forever)
    bnd_above = [];
end    

% RHS column
b_ent = [];
b_rows = [];
if ~isempty(bnd_above)
    b_ent = [b_ent;
        mu./4./pi./sqrt(r.^2 + (bnd_above - zp).^2);
        mu./4./pi.*(zp - bnd_above)./(sqrt(r.^2 + (bnd_above -zp).^2).^3)];
    b_rows = [b_rows, bnd_a_ind+(1:(2*rN))];
end
if ~isempty(bnd_below)
    b_ent = [b_ent;
        -mu./4./pi./sqrt(r.^2 + (bnd_below - zp).^2);
        -mu./4./pi.*(zp - bnd_below)./(sqrt(r.^2 + (bnd_below -zp).^2).^3)];
    b_rows = [b_rows, bnd_b_ind+(1:(2*rN))];
end

b = sparse(b_rows,1,b_ent,length(A),1);

% Back substitute
x = A \ b;
x = reshape(x,rN,[]);

% Extract the answer
ind1 = find(z_val<=z_range(1),1,'last');
if numel(z_range) == 2
    ind2 = find(z_val<=z_range(2),1,'last');
else
    ind2 = ind1;
end

% Output
xout = x(:,ind1:ind2);
rout = r;
zout = z_val(ind1:ind2);

    
    function [Ar r] = blk_ar()
    %----------------------------------------------------------------------
    % BLK_AR Generates one block of rN x rN matrix for the PDE
    %         d2/dr2 + 1/r d/dr A = 0 
    % with non-uniform, polynomially spaced grid.
    %
    % The function implements a Neumann condition at r=0, and a Dirichlet
    % condition at r=rbnd
    %
    % Non-uniform FDM techniques were used according to 
    %   Eugenia Kálnay de Rivas, "On the use of nonuniform grids in 
    %   finite-difference equations", Journal of Computational Physics, 
    %   Volume 10, Issue 2, October 1972, Pages 202-210 
    %
    % Inputs:
    %      rN - the number of elements to consider (global variable)
    %      rbnd - the upper bound of elements (global variable)
    %      ord - the order of the polynomial relationship (global variable)
    %
    % Outputs: 
    %      Ar - the single rN x rN FDM A matrix for the r variation.
    %      r  - the coordinates of each r point
    %----------------------------------------------------------------------
        xi = linspace(0,1,rN).'; % Define a uniformly spaced variable xi.
        dx = xi(2)-xi(1); % Extract xi spacing.
        r = (xi.^ord + 0.5./(rN-1).^ord) * rbnd; % Shift r to enforce Neumann
        dr = r(2)-r(1); % Extract r spacing.

        % 1 x rN row matrices of the entries involved.
        % r_fwd and r_rev to deal with the second derivative
        % r_dif to deal with the first derivative
        r_fwd = 1./(ord*(xi+dx/2).^(ord-1).*(ord*xi.^(ord-1)).*dx^2.*rbnd^2);
        r_rev = 1./(ord*(xi-dx/2).^(ord-1).*(ord*xi.^(ord-1)).*dx^2.*rbnd^2);
        r_dif = 1./(2*rbnd*dx .* r .* (ord*xi.^(ord-1)));

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
    end

        
    function A = blk_asm(Ar, zN, dz, psn)
    %--------------------------------------------------------------------------
    % BLK_ASM Assembles the full block of rNzN x rNzN matrix, with an 
    % optional poisson term. Implements the PDE
    %          Ar + d2/dz2 A + psn A = 0
    %
    % Implements the Dirichlet condition on both extremes of the z axis.
    % 
    % Inputs:
    %      Ar - the r-variational rNxrN FDM A matrix
    %      zN - the number of z elements to consider
    %      dz - the uniform separation distance between z elements
    %      psn - the poisson term (optional)
    %
    % Outputs: 
    %      A - the full rNzN x rNzN FDM A matrix for the two dimensional 
    %           block.
    %----------------------------------------------------------------------
        B = kron([1 -2 1]./dz^2,ones(zN,1));
        Az = spdiags(B,[-1 0 1],zN,zN);    
    
        A = kron(speye(zN),Ar)+kron(Az,speye(rN));
        if nargin == 4 && psn ~= 0
            A = A + psn*speye(rN*zN);
        end

    end

    function A = blk_bnd(A1,A2,fac)
    %--------------------------------------------------------------------------
    % BLK_BND Implements the continuity and gradient continuity boundary 
    % conditions:
    %           A2(1)               -  A1(end)                   = 0
    %          (A2(2)-A2(1))*fac(2) - (A1(end)-A1(end-1))*fac(1) = 0
    %
    % Can be used to implement a material boundary for both electrostatic
    % and magnetostatic problems
    %
    % Assumes:
    %      - Ar is common to both blocks and there is no mismatch of points
    %        along the r line. rN r elements are used.
    %      - The z axis is the second degree of freedom, and is uniformly
    %        spaced
    %      - Both matrices are square.
    % 
    % Inputs:
    %      A1 - the first matrix, which is being connected at its end
    %      A2 - the second matrix, being connected at its beginning
    %      fac - the factors used to resolve the gradient continuity
    %            condition
    %      rN - the number of r elements (global variable)
    % Outputs: 
    %      A - the full FDM A matrix for the two dimensional 
    %           block.
    %----------------------------------------------------------------------
        % Retrieve the number of z elements
        zN1 = length(A1)/rN; % number of z elements in A1
        zN2 = length(A2)/rN; % number of z elements in A2
    
        % Retrieve the d2/dz2 factors
        z_fac1 = A1(end, end-rN);
        z_fac2 = A2(1, 1+rN);
        
        % Make the boundary matrices
        Arows = zN1+[0,1,1,2,2,2,2,3];
        Acols = zN1+[1,1,2,0,1,2,3,2];
        Aents = [z_fac1, ...
                -1,1, ...
                [1,-1]*fac(1) , [-1,1]*fac(2), ...
                z_fac2];
        Abnd = sparse(Arows,Acols,Aents,zN1+zN2+2,zN1+zN2+2);
        Abnd = kron(Abnd, speye(rN));

        % Combine
        A = blkdiag(A1, sparse(rN*2,rN*2), A2);
        A = A + Abnd;
    end

end


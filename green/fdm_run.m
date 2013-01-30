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

% Error check before doing all the work!
assert(all(z_range<=bnds(coil_layer+1)), 'z goes out of upper bound!')
assert(all(z_range>=bnds(coil_layer)), 'z goes out of lower bound!')

% Make A block for each layer region
Alay = cell(layerN,1);
psn = -1j*w*mu.*sig.*mu_r; % poisson factors

% Make a vector of z values and boundaries
z_cell = cell(1,layerN);
bnd_ind = zeros(layerN,2); % index of the upper and lower boundaries
counter = 0;

[Ar r] = blk_ar; % radial blocks
for ii = 1:layerN 
    z_cell{ii} = linspace(bnds(ii),bnds(ii+1),zN(ii));
    bnd_ind(ii,1) = counter+1; % first index in this layer
    bnd_ind(ii,2) = counter+zN(ii); % last index in this layer
    counter = counter + zN(ii);
    
    % Create domain blocks
    Alay{ii} = blk_asm(Ar,z_cell{ii},psn(ii));
end

% Enforce boundary conditions
fac = 1./mu_r; % boundary factors
A = Alay{1};
for ii = 2:layerN
    A = blk_bnd(A,Alay{ii},z_cell{ii-1},z_cell{ii},fac(ii-1),fac(ii));
end
A = blk_dirichlet(A);

% DEBUG, make sure that z_val and A match or else all index extraction are
% wrong
assert((length(A)/rN) == counter);

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
        1./sqrt(r.^2 + (bnd_above - zp).^2);
        1.*(zp - bnd_above)./(sqrt(r.^2 + (bnd_above -zp).^2).^3)];
    b_rows = [b_rows, bnd_a_ind+(1:(2*rN))];
end
if ~isempty(bnd_below)
    b_ent = [b_ent;
        -1./sqrt(r.^2 + (bnd_below - zp).^2);
        -1.*(zp - bnd_below)./(sqrt(r.^2 + (bnd_below -zp).^2).^3)];
    b_rows = [b_rows, bnd_b_ind+(1:(2*rN))];
end

b = sparse(b_rows,1,b_ent,length(A),1);

% Back substitute
x = A \ b;
x = reshape(x,rN,[]);

% Extract the results in the coil layer
xout = x(:,bnd_ind(coil_layer,1):bnd_ind(coil_layer,2));
xout = xout*1e-7; %normalize to mu/4pi

% Interpolate
rout = r;
xout = interp1(z_cell{coil_layer},xout.',z_range,'cubic');
xout = xout.';
zout = z_range;

    
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
        xi = linspace(0,1,rN).'; % Define a uniformly space variable xi.
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
        B = [r_rev  -  r_dif, ...
             -r_rev -  r_fwd, ...
             r_fwd +  r_dif];

        % Form (rN-1) x rN matrix without x=end+1 row to implicitly enforce Dirichlet.
        Ar = spdiags(B(2:end,:), ...
            [0 1 2], rN-1,rN);

        % Enfore Neuman condition at x=0. Construct the matrix without
        % that row, then ENFORCE: x(0) = x(1).
        bnd_factor = 1/dr^2+1/(2*dr*r(1));
        neu_bnd = sparse([1; 1],[1; 2],...
            [-1; 1]*bnd_factor,...
            1,rN);
        Ar = [neu_bnd;Ar];
    end

        
    function A = blk_asm(Ar, z, psn)
    %--------------------------------------------------------------------------
    % BLK_ASM Assembles the domain block of rNzN x rNzN matrix, with an 
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
        zN1 = length(z);
    
        % Calculate the coefficients
        B = zeros(zN1-2,3);
        for ij = 1:length(B)
            B(ij,:) = fdcoeffF(2,z(ij+1),z(ij:ij+2));
        end
        Az = spdiags(B,[0 1 2],zN1-2,zN1);    
        
        % Indentity matrices to be kroneckered
        Areye = spdiags(ones(zN1-2,1),1,zN1-2,zN1);
        Azeye = speye(rN);
        
        A = kron(Areye,Ar)+kron(Az,Azeye);
        if nargin == 3 && psn ~= 0
            A = A + psn*kron(Areye,Azeye);
        end

    end

    function A = blk_bnd(A1,A2,z1,z2,fac1,fac2)
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
        
        % Continuity boundary condition
        %             zN1 zN2
        % [0 0 ... 0  -1   1  0 ... 0 0]
        Abnd1 = sparse([1,1],zN1+[0,1],[-1,1],1,zN1+zN2);
        Abnd1 = kron(Abnd1, speye(rN));

        % Gradient continuity boundary condition
        %             zN1 zN2
        % [0 0 ... 1  -1  -1  1 ... 0 0]
        ORDER = 2; % Gradient prediction FDM polynomial order
        sel_z1 = z1(end-ORDER:end); 
        sel_z2 = z2(1:1+ORDER);
        
        Arows = ones(2*(ORDER+1),1);
        Acols = zN1+(-ORDER:ORDER+1);
        c1 = fdcoeffF(1,sel_z1(end),sel_z1); % coefficients for d/dz A1
        c2 = fdcoeffF(1,sel_z2(1),sel_z2); % coefficients for d/dz A2
        Aents = [-c1*fac1, c2*fac2];
        
        Abnd2 = sparse(Arows,Acols,Aents,1,zN1+zN2);
        Abnd2 = kron(Abnd2, speye(rN));
        
        % Combine
        A = [A1, sparse(rN*(zN1-2),rN*zN2);
            Abnd1;
            Abnd2;
            sparse(rN*(zN2-2), rN*zN1), A2];
    end
    
    function A = blk_dirichlet(A)
    %----------------------------------------------------------------------
    % BLK_DIRICHLET Caps the top and bottom row of A with Dirichlet
    % conditions
    %----------------------------------------------------------------------
        totzN = length(A)/rN;
        Abnd1 = spdiags(ones(rN,1),0,rN,length(A));
        Abnd2 = spdiags(ones(rN,1),(totzN-1)*rN,rN,length(A));
        A = [Abnd1; A; Abnd2];
    end
end


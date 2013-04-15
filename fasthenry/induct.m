function [ V ] = induct(varargin)
%INDUCT: calculates the PEEC matrix vector product [L]*I for a flat
% coil systems, using the pre-corrected FFT algorithm.
%   
%   Usage: 
%   1a. Initialize the solver grid and kernel by calling
%           induct(xy_bnd,z_bnd); 
%       where xy_bnd is the maximum absolute bound for the x and y axes and
%       z_bnd is the bound for the z axis. 
%       
%       Initialize to an unspecified grid size by calling
%           induct(O,P);
%       a grid is then constructed for approximately 1 grid cell per 3
%       elements.
%
%   2.  For each new set of filaments, call
%           V = induct(I,O,L,W,H); or
%           induct(I,O,L,W,H);
%       to calculate the vector product [L]*I for [L] defined by the vectors
%       O, L, W, and H.
%
%   3.  To do more MV products, for example for GMRES iterations, call
%           V = induct(I);
%
%   Usage for subtracted Green's functions:
%   1b. Initialize the 1/r component and a subtract component by calling
%           induct(xy_bnd,z_bnd,tph_obj); 
%       where tph_obj is the subtracted Green's function object.
%
%   1b. Initialize the subtracted component without 1/r by calling
%           induct(xy_bnd,z_bnd,tph_obj,true); 
%       where tph_obj is the subtracted Green's function object.

persistent p
persistent nodes
persistent qwei
persistent U
persistent do_precor

% Internal settings
GRID_ELE_F = [2^6 2^6 2^2];
GRID_ELE_S = [2^6 2^6 2^8];
INTERP_ORD = 1;
DIRECT_STEN = 3;


% Initialization routine for free-space elements
if (nargin == 2 || nargin == 3 ||nargin==4) && nargout == 0
    % Plus or minus bounds where elements are contained
    xy_bnd = varargin{1};
    z_bnd = varargin{2};
    
    % Extract grid parameters if they are specified. 
    if numel(xy_bnd) == 1 && (numel(z_bnd)<=2)
        assert(numel(xy_bnd) == 1,'Only one value of xy_bnd accepted.');
        bnds = [-1 1; -1 1]*xy_bnd;
        if numel(z_bnd) == 1
            bnds = vertcat(bnds,z_bnd*[-1 1]);
        elseif numel(z_bnd) == 2
            bnds = vertcat(bnds,z_bnd(:)');
        else
            error('Either one absolute value or two complementary values for z_bnd.');
        end
        
        % Decide between the two predefined grid sizes.
        if nargin == 2
            g_e = GRID_ELE_F;
        else
            g_e = GRID_ELE_S; 
        end
    else
        % Else, xy_bnd is the lower left corner of each box, and z_bnd is
        % the upper right corner.
        N = size(xy_bnd,1);
        assert(size(z_bnd,1) == N);
        assert((size(z_bnd,2) == 3) && (size(xy_bnd,2) == 3));

        % Find the bounding box of the elements
        min_d = min([xy_bnd;z_bnd],[],1);
        max_d = max([xy_bnd;z_bnd],[],1);
        dist_d = max_d-min_d;

        % Find the number of boxes to make
        num_box(1) = ((N/3) * dist_d(1)^2/(dist_d(2)*dist_d(3)))^(1/3);
        num_box(2) = ((N/3) * dist_d(2)^2/(dist_d(1)*dist_d(3)))^(1/3);
        num_box(3) = ((N/3) * dist_d(3)^2/(dist_d(2)*dist_d(1)))^(1/3);
        %g_e = pow2(max(nextpow2(num_box)-1,0))*3; % Minimum 3 points!
        g_e = GRID_ELE_F;
        
        box_d = dist_d./(g_e-2)*1.005; % pad the boundary and make the borders slightly larger
        bnds = [min_d-box_d; max_d+box_d]';
    end
    % Initialize PFFT object
    p = pfft(g_e,bnds,INTERP_ORD,DIRECT_STEN,0);
    
    % Free-space kernel
    if nargin ~= 4
        p.init_kernel(@(x) 1e-7*inv_r(x));
        do_precor = true;
    else
        do_precor = false;
    end
    
    % Subtracted kernel
    if nargin > 2
        % Process TPH object
        tph_obj = varargin{3};
        assert(isa(tph_obj,'tph'),'Must supply a tph object!');

        % Setup green's functions
        t_fun = @(in) tph_wrap(tph_obj,'T',in);
        h_fun = @(in) tph_wrap(tph_obj,'H',in);
        
        % Init the TPH kernel
        p.init_kernel(t_fun, h_fun, [0 0 1]);
    end
    
    % Initialize quadrature rule
    [nodes, qwei] = nwspgr('KPU', 3, 3);
    
    fprintf('Initialization complete for g_e=%s\n',...
        mat2str(g_e));
    
    
% Calculation routine
elseif nargin == 5 && nargout <= 1
    a=tic;
    
    I = varargin{1};
    O = varargin{2};
    L = varargin{3};
    W = varargin{4};
    H = varargin{5};
    
    % Process Geometry
    % Setup the vector MV components
    L_mag = sqrt(sum(L.^2,2));
    U = bsxfun(@rdivide,L,L_mag); % unit vectors

    % Setup the centroids
    cen = O + L/2 + W/2 + H/2;

    % Setup quadrature points using sparse grid quadrature
    Nq = size(nodes,1);
    qpnts = cell(Nq,1);
    for ii = 1:Nq
        qpnts{ii} = O + bsxfun(@times,L,nodes(ii,:)) + bsxfun(@times,W,nodes(ii,:)) ...
                    + bsxfun(@times,H,nodes(ii,:));
    end
    
    % Project the geometry
    p.init_geometry(cen,qpnts,qwei,L_mag);
    fprintf('Pfft projection time: %g\n',toc(a));
    
    clear cen qpnts qwei L_mag
    if do_precor
        a = tic;
        
        % Setup the index-based calcp function
        getidx = @(mat,from,to) [mat(from,:) mat(to,:)]';
        this_calcp = @(from,to) mutual(getidx(O,from,to), getidx(L,from,to),...
        getidx(W,from,to), getidx(H,from,to));
        
        p.init_precorrect(this_calcp,true);
        fprintf('Pfft precorrect time: %g\n',toc(a));
    end
    
    % Deliver answer
    if nargout > 0
        V = p.fastmv(I,U);
    end
    
% Single MV product
elseif nargin == 1 && nargout == 1 && any(any(U))
    V = p.fastmv(varargin{1},U);
    
else
    help induct_free;
end

end


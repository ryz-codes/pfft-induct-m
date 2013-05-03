function [x,y,z,E,B] = fastfield(g_e,bnds,geom,Ib,sub)
%FASTFIELD Uses pFFT to project the current out to a grid, and use pFFT to
%perform convolutions and spit out the fields.
% rad is the radius of the path. 
% sub dictates whether the add the 1e-7/r component (default is yes)
O = geom.O; L = geom.L; W = geom.W; H = geom.H;

% Check the problem type.
if (iscell(O) && iscell(L) && iscell(W) && iscell(H))
    if (iscell(O{1}) && iscell(L{1}) && iscell(W{1}) && iscell(H{1}))
        TYPE_ = 2;
    else
        TYPE_ = 1;
    end
else
    TYPE_ = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch TYPE_ % 0 no terminal, 1 single terminal, 2 multi terminal
%--------------------------------------------------------------------------    
    case 0 % No terminal
        error('No-terminal extraction not implemented');
    case 1 % 1 terminal
        % Expand out the branches
        O = vertcat(O{:}); L = vertcat(L{:}); W = vertcat(W{:}); H = vertcat(H{:});
        % Terminal matrix
        S = 1; %one terminal
    case 2 % multi-terminal
        Nterms = numel(O);
        for ij = 1:Nterms
            % Expand out the branches
            O{ij} = vertcat(O{ij}{:}); L{ij} = vertcat(L{ij}{:}); 
            W{ij} = vertcat(W{ij}{:}); H{ij} = vertcat(H{ij}{:});
        end
        % Expand out the terminals
        O = vertcat(O{:}); L = vertcat(L{:}); 
        W = vertcat(W{:}); H = vertcat(H{:});
end

% Initiate PFFT
p = pfft(g_e,bnds,1,2,0);

if isfield(geom,'g')
    tph_obj = geom.g;
    assert(isa(tph_obj,'tph'),'Must supply a tph object!');
    
    % clear out the kernels.
    p.init_kernel; 
    
    % Setup green's functions
    if nargin >= 5 && sub == true
        t_fun = @(in) tph_wrap(tph_obj,'T',in); % Just Toeplitz
        fprintf('1/r not added\n');
    else
        t_fun = @(in) tph_wrap(tph_obj,'TT',in); % Toeplitz plus 1/r
    end
    h_fun = @(in) tph_wrap(tph_obj,'H',in);

    % Init the new TPH kernel
    p.init_kernel(t_fun, h_fun, [0 0 1]);
else
    p.init_kernel(@(x) 1e-7*inv_r(x));
end

% Perform projection
a=tic;

% Process Geometry
% Setup the vector MV components
L_mag = sqrt(sum(L.^2,2));
U = bsxfun(@rdivide,L,L_mag); % unit vectors

% Setup the centroids
cen = O + L/2 + W/2 + H/2;

% Setup quadrature points using sparse grid quadrature
[nodes, qwei] = nwspgr('KPU', 3, 3);
Nq = size(nodes,1);
qpnts = cell(Nq,1);
for ii = 1:Nq
    qpnts{ii} = O + bsxfun(@times,L,nodes(ii,:)) + bsxfun(@times,W,nodes(ii,:)) ...
                + bsxfun(@times,H,nodes(ii,:));
end

% Project the geometry
p.init_geometry(cen,qpnts,qwei,L_mag);
fprintf('Pfft projection time: %g\n',toc(a));
    
[~, out_proj] = p.fastmv(Ib,U); % Electric field in each direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ALL THE FIELDS
% x,y,z
x = linspace(min(bnds(1,:)),max(bnds(1,:)),g_e(1));
y = linspace(min(bnds(2,:)),max(bnds(2,:)),g_e(2));
z = linspace(min(bnds(3,:)),max(bnds(3,:)),g_e(3));

% Electric field
E = cell(1,4);
E{4} = zeros(g_e);
for ii = 1:3
    if ii <= size(out_proj,2)
        E{ii} = reshape(out_proj(:,ii),g_e);
        % convert data from ndgrid to meshgrid
        E{ii} = permute(E{ii},[2 1 3]);
        E{4} = E{4} + abs(E{ii}).^2;
    else
        E{ii} = zeros(g_e);
    end
end
E{4} = sqrt(E{4}); % Electric field strength

% Magnetic field is curl of E
B = cell(1,4); B{4} = 0;
[gx,gy,gz] = meshgrid(x,y,z);
[B{1},B{2},B{3}] = curl(gx,gy,gz,E{1},E{2},E{3});
for ii = 1:3
    B{4} = B{4} + abs(B{ii}).^2;
end
B{4} = sqrt(B{4});

% stream line plot
num_lines = 5;
x_start = linspace(min(x),max(x),num_lines);
y_start = linspace(min(y),max(y),num_lines);
z_start = linspace(min(z),max(z),num_lines);
[sx,sy,sz] = meshgrid(x_start,y_start,z_start);
clf;
h1 = streamline(gx,gy,gz,E{1},E{2},E{3},sx,sy,sz);
h2 = streamline(gx,gy,gz,B{1},B{2},B{3},sx,sy,sz);
set(h1,'Color','blue'); set(h2,'Color','red');

end


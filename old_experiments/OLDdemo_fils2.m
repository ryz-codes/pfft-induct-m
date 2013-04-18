

clear
load ferrite_iron.mat;

RAD = 0.1;
GRID_ELE = [2^6 2^6 2^4];
OUTER_BOUNDS = [[-1 1; -1 1]*RAD*1.2; [-0.01 -0.001]]; % Plus or minus bounds where elements are contained
INTERP_ORD = 2;
DIRECT_STEN = 4;

% Setup green's functions
t_fun = @(in) tph_wrap(g,'T',in);
h_fun = @(in) tph_wrap(g,'H',in);

p = pfft(GRID_ELE,OUTER_BOUNDS,INTERP_ORD,DIRECT_STEN);
p.init_kernel(@(x) 1e-7*inv_r(x));
p.init_kernel(t_fun, h_fun, [0 0 1]);
%% Setup Geometry

fprintf('\n\n');
NUM_ELE = 300; % Number of elements simulated
WID = 0.001;
HEI = 0.001;

% Generate a circle of filaments
[O L W H] = genCircFils( RAD, -0.005, NUM_ELE, WID, HEI );
I = ones(NUM_ELE,1);

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

% Setup the index-based calcp function
getidx = @(mat,from,to) [mat(from,:) mat(to,:)];
this_calcp = @(from,to) mutual(getidx(O,from,to)',... 
    getidx(L,from,to)', getidx(W,from,to)', getidx(H,from,to)');

% Run pfft

p.init_geometry(cen,qpnts,qwei,L_mag);
clear cen qpnts qwei L_mag
p.init_precorrect(this_calcp,true);
Vpfft = p.fastmv(I,U);
a = sum(Vpfft);
fprintf('Pfft solution: %g\n',a);

% Run direct solve
if NUM_ELE < 800
    Vcalcp = calcp(I,O,L,W,H);
    b = sum(Vcalcp);
    fprintf('Calcp solution: %g\n',b);
    fprintf('Pfft error: %g\n',abs(a-b)/abs(b));
end
NUM_ELE = 100000; % Number of elements simulated
WID = 5e-4;
HEI = 5e-4;
RAD = linspace(2.5e-2,10.5e-2,23);

% Generate a circle of filaments
[O L W H] = genCircFils( RAD, 0, ceil(NUM_ELE/numel(RAD)), WID, HEI );

xy_bnd = O; z_bnd = O + L + W + H;

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
g_e = pow2(max(nextpow2(num_box)-2,0))*3; % Minimum 3 points!

box_d = dist_d./(g_e-2);
bnds = [min_d-box_d; max_d+box_d]';

disp(g_e)
disp(bnds)
function [fils_struct, Ntotfils] = genCoilLitzFils(refine, N_fil, pxrad,z_off)
%GENCOILFILS Generates filaments used to model the coil used in our 
% experiments
%
% The resultant coil is centered at x=0 and y=0.
%
% Measurement errors are about +- 0.5mm
% Relative max error about 0.0307
% Average meas error about 0.02
% DEFAULTS
% Correct fill factor is 16* pi * (wdiam/2)^2/wei/hei.
if nargin < 1
    refine = [1 2];
end
if nargin < 2
    N_fil = 1000;
end
if nargin < 3
    pxrad = 4;
end
if nargin <4
    z_off=0;
end


% Lisa's measurements
% first set of turns
a1(1)= 1+15/16-14/32; %inch
a1(2)= 1+15/16;
a1 = a1*25.4;
N1 = 7;

% second set of turns
a2(1) = 3.25-18/32;
a2(2) = 3.25;
a2 = a2*25.4;
N2 = 9;

% third set of turns
a3(1) = 4.5-22/32;
a3(2) = 4.5;
a3 = a3*25.4;
N3 = 12;

% the wires themselves
wid = (0.063*25.4)*1e-3; %[m] 0.063 in +- 0.002in
hei = 7.493e-3; %[m] 0.295 in +- 0.005in
z_off = z_off - hei/2;

% turn into sampled points
rho = [a1,a2,a3];
theta(1) = 0;
theta(2) = N1*2*pi;
theta(3) = theta(2) + pi/2; 
theta(4) = theta(3) + N2*2*pi - pi/4;
theta(5) = theta(4) + pi/4;
theta(6) = theta(5) + N3*2*pi - pi/2;

% Interpolate between them to get a function of rho_c wrt to theta
% Do this to maintain a constant filament length
lens = cumtrapz(theta,rho);
lens_c = linspace(lens(1),lens(end),N_fil+1);
theta_c = interp1(lens,theta,lens_c,'linear');
rho_c = interp1(lens,rho,lens_c,'linear')*1e-3;
    
coil_path(:,1) = rho_c.*cos(theta_c);
coil_path(:,2) = rho_c.*sin(theta_c);
coil_path(:,3) = z_off*ones(size(theta_c));

fprintf('Average filament length: %g\n',mean(sqrt(sum(diff(coil_path,1).^2,2))));

%---
% NOW Generate the filaments and divide up the path
%---
filcoil = genFilsFromPath(coil_path,wid,hei);
[O, L, W, H] = refineCrossSec(refine(1),refine(2),filcoil);

% Begin to twist the litz
num_litz = prod(refine);
fils = permute(reshape([O L W H]',12,N_fil,num_litz),[2 1 3]);
% matrix ordering: filaments, dimensions, strand

% Shift the ownership of the bundles.
for ij = 1:N_fil
    shft = randperm(num_litz);
    %shft = mod((1:N_fil)-1,num_litz) + 1;

    % Randomize ownership of the bundle. 
    fils(ij,:,1:num_litz) = fils(ij,:,shft);
end

Ntotfils = numel(fils(:,1,:));

% For each litz strand, refine and store as one segment
O = cell(num_litz,1); L = cell(num_litz,1);
W = cell(num_litz,1); H = cell(num_litz,1);
for ii = 1:num_litz
    O{ii} = cell(N_fil,1);
    L{ii} = cell(N_fil,1);
    W{ii} = cell(N_fil,1);
    H{ii} = cell(N_fil,1);
    for ij = 1:N_fil
        [O{ii}{ij},L{ii}{ij},W{ii}{ij},H{ii}{ij}]...
            = fil2fils(squeeze(fils(ij,1:3,ii)),...
                       squeeze(fils(ij,4:6,ii)), ...
                       squeeze(fils(ij,7:9,ii)),...
                       squeeze(fils(ij,10:12,ii)),pxrad);
        if nargin == 0
          showFils(O{ii}{ij},L{ii}{ij},W{ii}{ij},H{ii}{ij});
          pause
        end
    end
end

Ntotfils = Ntotfils * size(O{1}{1},1);

fils_struct.O = O;
fils_struct.L = L;
fils_struct.W = W;
fils_struct.H = H;
% if nargin == 0
%     %showFils(vertcat(O{1}{:}),vertcat(W{1}{:}),vertcat(L{1}{:}),vertcat(H{1}{:}))
%     view([180-45 35]);
% end
end


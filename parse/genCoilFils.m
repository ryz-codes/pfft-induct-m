function [fils] = genCoilFils(REFINE_LEV, N_fil, z_off)
%GENCOILFILS Generates filaments used to model the coil used in our 
% experiments
%
% The resultant coil is centered at x=0 and y=0. The top of the coil is
% placed at z=z_off.
%
% Measurement errors are about +- 0.5mm
% Relative max error about 0.0307
% Average meas error about 0.02
% DEFAULTS
if nargin < 1
    REFINE_LEV = 1;
end
if nargin < 2
    N_fil = 600;
end
if nargin < 3
    z_off = 0;
end
BRIDGELEN = 0.01; %m

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

% % My meas
% % first set
% a1(2) = 3.85*25.4/2;
% a1(1) = a1(end) - 10;
% N1 = 7;
% 
% % second set
% a2(1) = a1(end) +20;
% a2(2) = a2(1) + 0.55*25.4;
% N2 = 9;
% 
% % third set
% a3(1) = a2(end) + 11;
% a3(2) =a3(1)+0.74*25.4;
% N3 = 12;
% 
% the wires themselves
wid = (0.063*25.4)*1e-3; %[m] 0.063 in +- 0.002in
hei = 7.493e-3; %[m] 0.295 in +- 0.005in

% the board
t_b = 0.017; %[m]

%z_off = z_off-hei/2;

% turn into sampled points
rho = [a1,a2,a3];
theta(1) = 0;
theta(2) = N1*2*pi;
theta(3) = theta(2) + pi/2; 
theta(4) = theta(3) + N2*2*pi - pi/4;
theta(5) = theta(4) + pi/4;
theta(6) = theta(5) + N3*2*pi - pi/2;

% interpolate between them!
theta_c = linspace(theta(end),theta(1),N_fil);
rho_c = interp1(theta,rho,theta_c,'linear')*1e-3; %[m]

coil_path(:,1) = rho_c.*cos(theta_c);
coil_path(:,2) = rho_c.*sin(theta_c);
coil_path(:,3) = zeros(size(theta_c));

%------------------
% BRIDGES
%------------------
% make into path and add bridges
B_FILLEN = rho_c(1) * (theta_c(1)-theta_c(2));

% outer path 
lil_b = linspace(rho_c(1)+BRIDGELEN,rho_c(1),ceil(1+BRIDGELEN/B_FILLEN));
lil_b = lil_b(1:end-1); %overlapping point
coil_path = [[lil_b(:), zeros(length(lil_b),2)]; coil_path];

% inner path 
lil_b = linspace(rho_c(end),0,ceil(1+rho_c(end)./B_FILLEN));
lil_b = lil_b(2:end); %overlapping point
coil_path = [coil_path; [lil_b(:), zeros(length(lil_b),2)]];

% bridge down
b = linspace(0,-t_b-hei,ceil(1+abs((-t_b-hei)./B_FILLEN)));
bx = linspace(0,1e-3,length(b));
bridge_path = [bx(:),zeros(length(b),1),b(:)];

% bridge back
b = linspace(1e-3,rho_c(1)+BRIDGELEN,ceil(1+(rho_c(1)+BRIDGELEN)./B_FILLEN));
b = b(2:end);
bridge_path = [bridge_path;[b(:),zeros(length(b),1),repmat(-t_b-hei,length(b),1)]];

%---
% NOW OFFSET IN Z
%---
coil_path(:,3) = coil_path(:,3) + z_off;
bridge_path(:,3) = bridge_path(:,3) + z_off;

filcoil = genFilsFromPath(coil_path,wid,hei);
filcoil = refineCrossSec(1,REFINE_LEV,filcoil);

filbridge = genFilsFromPath(bridge_path,hei,wid);
filbridge = refineCrossSec(REFINE_LEV,1,filbridge);

fils = combineFils(filcoil,filbridge);

if nargin == 0
    showFils(fils)
    view([180-45 35]);
end
end


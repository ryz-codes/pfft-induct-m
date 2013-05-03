function [coil_path] = genCoilPath(N_fil)
%GENCOILFILS Generates filaments used to model the coil used in our 
% experiments
%
% The resultant coil is centered at x=0 and y=0.
%
% Measurement errors are about +- 0.5mm
% Relative max error about 0.0307
% Average meas error about 0.02
% DEFAULTS
if nargin < 1
    N_fil = 1000;
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
coil_path(:,3) = zeros(size(theta_c));

%---
% NOW Generate the filaments and divide up the path
%---
coil_path;


if nargin == 0
    showPath(coil_path)
end
end


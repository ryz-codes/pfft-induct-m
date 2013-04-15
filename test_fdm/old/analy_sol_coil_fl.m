function [A] = analy_sol_coil_fl(R,z,L)
% This code compares the data extracted from the FDM simulation with the
% analytical solution commonly found in literature. The analytical solution
% is generated using the fast hankel transform
% 
persistent k1
persistent r1
persistent I
%% Prepare Kernel
if isempty(I)
    [~,k1,r1,I]=fht(@(x) 1,R*10,2*pi*100/R,1);
    disp('Recalculating kernel for FL');
end

%% Preparation
c = L.coil_layer;

% plate thicknesses
c_bnd = L.bnds(c);
t(1) = c_bnd - L.bnds(c-1); % bottom plate thickness
t(2) = L.bnds(c+2) - L.bnds(c+1); % top plate thickness
s = L.bnds(c+1) - c_bnd; % gap thickness

mu=1e-7*4*pi;
w = L.w;
mu_r(1) = L.mu_r(c-1);
sig(1) = L.sig(c-1);
mu_r(2) = L.mu_r(c+1);
sig(2) = L.sig(c+1);

% Setup up components
eta = @(k,w,id) sqrt(k.^2+1j*w*mu*mu_r(id)*sig(id));
phi = @(k,w,id) (mu_r(id)*k-eta(k,w,id))./(mu_r(id)*k+eta(k,w,id)); % Halfspace Green's fn
lambda = @(k,w,id) phi(k,w,id) .* (1-exp(-2*t(id)*eta(k,w,id))) ./ ...
    (1-phi(k,w,id).^2.*exp(-2*t(id)*eta(k,w,id))); % Single layer Green's fn

% Get parameters
d1 = z-c_bnd;
d2 = z-c_bnd;
d1p = s-d1;
d2p = s-d2;

% Set up Green's function
f = @(k,w) (lambda(k,w,1).*exp(-k.*(d1+d2)) + lambda(k,w,2).*exp(-k.*(d1p+d2p))) ...
    ./ (1 - lambda(k,w,1).*lambda(k,w,2).*exp(-2*k*s)); % Hankel part
g = @(k,w) (2.*lambda(k,w,1).*lambda(k,w,2).*exp(-2*k*s).*cosh(k.*(d2-d1))) ...
    ./ (1 - lambda(k,w,1).*lambda(k,w,2).*exp(-2*k*s)); % Toeplitz part

% Setup analytical Green's function
Ar_ker = @(k,z,w) mu * pi * R^2 ./ k .* ((f(k,w)+g(k,w))) .* besselj(1,k*R); 

% Generate results at z=0 and at z=0.05 for d = 0.01.
% A1 is at z=d, A2 is at z=0
% mu, mu_r, d, w should all be included in FDM_data.mat
A = ifht(Ar_ker(k1,z,w),k1,r1,I);
A = interp1(r1,A,R);
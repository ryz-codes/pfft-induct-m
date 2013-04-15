function [A r] = analy_sol_hs(zh,L)
% This code compares the data extracted from the FDM simulation with the
% analytical solution commonly found in literature. The analytical solution
% is generated using the fast hankel transform

persistent k1
persistent r1
persistent I
%% Prepare Kernel
if isempty(I)
    [~,k1,r1,I]=fht(@(x) 1,5,10000,0);
end

zh = abs(zh);
mu=1e-7*4*pi;
w = L.w;
mu_r = max(L.mu_r);
sig = max(L.sig);

% Setup analytical Green's function
eta = @(k,w) sqrt(k.^2+1j*w*mu*mu_r*sig);
phi = @(k,w) (mu_r*k-eta(k,w))./(mu_r*k+eta(k,w));
Gr_ker = @(k,z,w) 1e-7 ./ k .* (phi(k,w).* exp(-k.*(z)));
            
% Generate results at z=0 and at z=0.05 for d = 0.01.
% A1 is at z=d, A2 is at z=0
% mu, mu_r, d, w should all be included in FDM_data.mat
A = ifht(Gr_ker(k1,zh,w),k1,r1,I);
r = r1;
function [A] = analy_coils_plate(R,z,L)
% This code compares the data extracted from the FDM simulation with the
% analytical solution commonly found in literature. The analytical solution
% is generated using the fast hankel transform
% 
persistent k1
persistent r1
persistent I
% %% Prepare Kernel
% if isempty(I)
%     addpath([pwd '/hankel']) % Include the package
%     [~,k1,r1,I]=fht(@(x) 1,5,10000,0);
% end

%% Preparation
z = 2*abs(z);
mu=1e-7*4*pi;
w = L.w;
mu_r = max(L.mu_r);
sig = max(L.sig);
plate_lay = find(L.sig == sig);
t = L.bnds(plate_lay+1) - L.bnds(plate_lay);

% Setup analytical Green's function
eta = @(k,w) sqrt(k.^2+1j*w*mu*mu_r*sig);
phi = @(k,w) (mu_r*k-eta(k,w))./(mu_r*k+eta(k,w));
lambda = @(k,w) phi(k,w) .* (1-exp(-2*t*eta(k,w))) ./ ...
    (1-phi(k,w).^2.*exp(-2*t*eta(k,w))); % Single layer Green's fn
            
% Generate the Kernals
if isempty(I)
    [~,k1,r1,I]=fht(@(x) 1,max(R)*10,2*pi*10/min(R),1);
end

% Superposition the contribution of the rings
%Ar_ker_num = 0;
Ar_ker_num = besselj(1,bsxfun(@times,k1(:).',R(:)));
Ar_ker_num = sum(bsxfun(@times,Ar_ker_num,R(:)),1);
Ar_ker_num = (mu * pi ./ k1 .* (lambda(k1,w).* exp(-k1.*(z)))) .* Ar_ker_num;

% Generate results at z=0 and at z=0.05 for d = 0.01.
% A1 is at z=d, A2 is at z=0
% mu, mu_r, d, w should all be included in FDM_data.mat
A = r1.*ifht(Ar_ker_num,k1,r1,I);
A = interp1(r1,A,R);
A = sum(A);
function [rel_err, abs_err] = test_error(z, zp, L)
%TEST_ERROR Tests the FDM error for the half-space and 5-plate
%configurations against analytical solutions
%

if nargin ==0
z = 1e-3;
zp = 5e-3;
L = defaultL(2);
end

%% Pick the right layer configuration
% Calculate using the FDM
fprintf('Running FDM...\n');
[A1 r] = fdm_run(z,zp,L);

%% Prepare Kernel
mu=1e-7*4*pi;
w = L.w;

% Pick the correct analytical expression for the Green's funciton of such a
% layered geometry.
switch L.layerN
    case 2 % HALFSPACE
        fprintf('Recalculate using Half-Space Analytical formula.\n');
        
        zh = abs(z+zp);
        mu_r = max(L.mu_r);
        sig = max(L.sig);

        % Setup analytical Green's functions
        eta = @(k,w) sqrt(k.^2+1j*w*mu*mu_r*sig);
        phi = @(k,w) (mu_r*k-eta(k,w))./(mu_r*k+eta(k,w));
        integ = @(k) 1 ./ k .* (phi(k,w).* exp(-k.*(zh)));
        
    case 5 % SANDWICH
        fprintf('Recalculate using Sandwich Analytical formula.\n');
        
        % coil layer
        c = L.coil_layer;

        % plate thicknesses
        c_bnd = L.bnds(c);
        t(1) = c_bnd - L.bnds(c-1); % bottom plate thickness
        t(2) = L.bnds(c+2) - L.bnds(c+1); % top plate thickness
        s = L.bnds(c+1) - c_bnd; % gap thickness

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
        d1 = zp-c_bnd;
        d2 = z-c_bnd;
        d1p = s-d1;
        d2p = s-d2;

        % Set up Green's function
        f = @(k,w) (lambda(k,w,1).*exp(-k.*(d1+d2)) + lambda(k,w,2).*exp(-k.*(d1p+d2p))) ...
            ./ (1 - lambda(k,w,1).*lambda(k,w,2).*exp(-2*k*s)); % Hankel part
        g = @(k,w) (2.*lambda(k,w,1).*lambda(k,w,2).*exp(-2*k*s).*cosh(k.*(d2-d1))) ...
            ./ (1 - lambda(k,w,1).*lambda(k,w,2).*exp(-2*k*s)); % Toeplitz part

        % Generate results at z=0 and at z=0.05 for d = 0.01.
        % A1 is at z=d, A2 is at z=0
        % mu, mu_r, d, w should all be included in FDM_data.mat
        integ = @(k) 1 ./ k .* (g(k,w) + f(k,w)); 
        
    otherwise % should never happen!
        error('Number of layers provided: %d. This is not supported.',L.layerN);
end
    
%%
% Inverse hankel transform using quadrature.
K = 1e4; % Maximum quadrature frequency.
tic
[A2, err] = quadht(integ,K,r,0,5e4); % use an intense quadrature to evaluate
fprintf('---QUADRATURE---\nQuadrature time:\t%g seconds\n',toc);
fprintf('Quadrature rel err:\t%1.1e\n',err/norm(A2,inf));
fprintf('Rel Err due to cutoff:\t%1.1e\n---\n',integ(K)/norm(A2,inf));
A2 = A2*1e-7; % renormalize to mu/4pi

% % The following plots the integ function
% k1 = K*logspace(-4,1);
% figure(1)
% loglog(k1,abs(real(integ(k1))),k1,abs(imag(integ(k1))));

%% graphical treatment
figure(2)
subplot(211)
loglog(r,abs(real(A1)),r,abs(real(A2)),r,abs(real(A1)-real(A2)));
subplot(212)
loglog(r,abs(imag(A1)),r,abs(imag(A2)),r,abs(imag(A1)-imag(A2)));


%% Error estimates
abs_err = norm(A2-A1,inf);
rel_err = abs_err/norm(A2,inf);
fprintf('---FDM--\nerr: %1.1e\trel err: %1.1e\n\n',abs_err,rel_err);
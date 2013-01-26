%function [A, r] = test_error(z, zp, def_num)

z = 1e-3;
zp = 5e-3;
def_num = 2;
%ANALY_GREEN Analytical solution for Green's function
%   Calculates the analytical solution for the Green's function of:
%     half-space (2 half-spaces and one boundary)
%     sandwich (3 intermediate layers, 4 boundaries)
%
%   To do:
%     plate (2 half-spaces, one intermediate "plate" layer, 2 boundaries)


%% Pick the right layer configuration
L = defaultL(def_num);
K = 1e4;
assert(def_num == 2 || def_num==5, 'Only defaultL(2) and defaultL(5) supported right now!');

% Calculate using the FDM
fprintf('Running FDM...\n');
[A1 r z_actual] = fdm_run(z,zp,L);
z = z_actual;

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
tic
[A2, err] = quadht(integ,K,r,0,5e4); % use an intense quadrature to evaluate
fprintf('---\nQuadrature time: %g seconds\n',toc);
fprintf('abs err: %1.1e\trel err: %1.1e\n---\n',err,err/norm(A2,inf));
A2 = A2*1e-7; % renormalize

%% graphical treatment
k1 = K*logspace(1e-4,1);
figure(1)
loglog(k1,abs(real(integ(k1))),k1,abs(imag(integ(k1))));

figure(2)
subplot(211)
loglog(r,abs(real(A1)),r,abs(real(A2)),r,abs(real(A1)-real(A2)));
subplot(212)
loglog(r,abs(imag(A1)),r,abs(imag(A2)),r,abs(imag(A1)-imag(A2)));

% Error estimates
fdmerr = norm(A2-A1,inf);
fprintf('FDM err: %1.1e\trel err: %1.1e\n\n',fdmerr,fdmerr/norm(A2,inf));

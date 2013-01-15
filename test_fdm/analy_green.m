function [A, r, rel_err, abs_err] = analy_green(z, zp, L)
%ANALY_GREEN Analytical solution for Green's function
%   Calculates the analytical solution for the Green's function of:
%     half-space (2 half-spaces and one boundary)
%     sandwich (3 intermediate layers, 4 boundaries)
%
%   To do:
%     plate (2 half-spaces, one intermediate "plate" layer, 2 boundaries)

persistent k1
persistent r1
persistent I1
%% Prepare Kernel
zh = abs(z+zp);
if isempty(I1)
    R = 0.1;
    K = 3e4;
    [~,k1,r1,I1]=fht([],R,K,0,4,6);
end

mu=1e-7*4*pi;
w = L.w;
switch L.layerN
    case 2 % HALFSPACE
        fprintf('Half-space configuration detected.\n');
        
        mu_r = max(L.mu_r);
        sig = max(L.sig);

        % Setup analytical Green's functions
        eta = @(k,w) sqrt(k.^2+1j*w*mu*mu_r*sig);
        phi = @(k,w) (mu_r*k-eta(k,w))./(mu_r*k+eta(k,w));
        Gr_ker = @(k) 1e-7 ./ k .* (phi(k,w).* exp(-k.*(zh)));
        inv_trans(Gr_ker);
        
    case 5 % SANDWICH
        fprintf('Sandwich configuration detected.\n');
        
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
        integ = @(k) 1e-7 ./ k .* (g(k,w) + f(k,w)); 
        inv_trans(integ);
        
    otherwise
        error('Number of layers provided: %d. This is not supported.',L.layerN);
       
end
    function inv_trans(fun)
    % Inverse hankel transform "fun" and derive estimates of the error.
        A = ifht(fun(k1),k1,r1,I1)+endcor(fun,k1(1),r1,0);
        r = r1;

        r2 = linspace(r1(1),r1(end),5); % pick 5 points to test.
        A2 = endcor(fun,k1(end)*1e1,r2,0,1e6); % use an intense quadrature to evaluate
        A2est = interp1(r1,A,r2);
        
        abs_err = max(abs(A2est-A2));
        rel_err = abs_err/max(abs(A2));
        if nargout < 3
            fprintf('rel_err: %1.1e\tabs_err: %1.1e\n\n',rel_err,abs_err);
        end
        
        % graphical treatment
        loglog(k1,abs(real(fun(k1))),k1,abs(imag(fun(k1))));
    end
end


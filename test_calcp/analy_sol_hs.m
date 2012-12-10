function [M] = analy_sol_hs(R,z,zp,w)
% This code gives the analytical solution for the subtracted
% self-inductance of a ring.

% Experimentally verified analytical solution is from Eqn (A23) of
%   Hurley, W.G.; Duffy, M.C.; , "Calculation of self and mutual impedances 
%   in planar magnetic structures," Magnetics, IEEE Transactions on , 
%   vol.31, no.4, pp.2416-2422, Jul 1995
%   doi: 10.1109/20.390151
%   URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=390151&isnumber=8841


%% Preparation
addpath([pwd '/hankel']) % Include the fht package

mu=1e-7*4*pi;
mu_r = 130;
sig = 8e6;
a = R;

% Setup analytical solution for a ring above an interface
eta = @(k,w) sqrt(k.^2+1j*w*mu*mu_r*sig);
phi = @(k,w) (mu_r*k-eta(k,w))./(mu_r*k+eta(k,w));
Gr_ker = @(k,z,zp,w) mu * pi * a ./ k .* besselj(1,k.*a).*...
    (exp(-k.*abs(z-zp))+phi(k,w).* exp(-k.*abs(z+zp)));
            
% Generate the Kernals
[temp,k1,r1,I]=fht(@(x) 1,1,10000,1);

%figure;loglog(k1,abs(real(Gr_ker(k1,z,zp,w))));

% Perform inverse hankel transform
A = ifht(Gr_ker(k1,z,zp,w),k1,r1,I);
A = A.*r1;

% extract mutual
M = interp1(r1,A,R);
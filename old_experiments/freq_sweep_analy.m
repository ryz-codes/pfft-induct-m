% FREQ_SWEEP_ANALY.M
% The following script produces a plot of the analytical solution versus
% numerically simulated for a fixed discretization and for varying frequencies. 
%
% Two plots are produced, one for L0, DeltaL (half-space) and 
% DeltaL(sandwich), and one for DeltaR and DeltaR.


clear all
addpath([pwd '/test_fdm'])
%% Setup the frequencies and geometries to be swept
%  Setup the Green's function if needed
w = logspace(2,6,10);
RECALC_G = 0; % Re-evaluate the green's functions
RECALC_Z = 0; % Re-evaluate the impedances

% Generate Filaments
NUM_ELE = 300; % Number of elements simulated
WID = 1e-6;
HEI = 1e-6;
D = 4e-3;
RAD = 0.1;%linspace(2.5e-2,10.5e-2,23);

% Generate circular filaments
path = genSpiralPath(RAD, -(D+HEI/2), NUM_ELE);
fils_hs = genFilsFromPath(path,WID,HEI);

path = genSpiralPath(RAD, (D+HEI/2), NUM_ELE);
fils_fl = genFilsFromPath(path,WID,HEI);

% Excitation current is uniform
I = ones(size(fils_hs{1},1),1);

% Setup the Green's functions
try 
    obj = load('g_swept','w');
catch
    obj = struct('w',[]);
end

if RECALC_G || ~isequal(obj.w,w) % Manual override or if w has changed
%-------------------------------------------------------------------------
% The following is the time-consuming evaluation of Green's functions.
% Perform once!
g_hs = sweep_tph(w,1); % Use defaultL(1).
g_fl = sweep_tph(w,5); % Use defaultL(5).
save('g_swept','g_hs','g_fl','w');
%-------------------------------------------------------------------------
else
    load('g_swept');
end

%% Calculate the Free-space component
% Setup the kernel
induct(max(RAD)*1.2,[-0.01 -0.002]);
Vfree = sum(induct(I,fils_hs{:}));
Vfreea = analySqRing(RAD,WID);
(Vfree-Vfreea)/Vfreea

%% Iterate through each frequency
if RECALC_Z
% Half-space
Vhs = zeros(length(w),1); % impedance calculated
Vhsa = zeros(length(w),1); % impedance analytical
% Five-layer
Vfl = zeros(length(w),1); % impedance calculated
Vfla = zeros(length(w),1); % impedance calculated

for ii = 1:length(w)
    %%% HALF-SPACE
    % pfft solution
    induct(max(RAD)*1.2,[-0.010 -0.002],g_hs{ii},true);
    Vhs(ii) = sum(induct(I,fils_hs{:}));
    
    % analytic solution
    thisL = defaultL(1); thisL.w = w(ii);
    Vhsa(ii) = analy_sol_coil(RAD,D,thisL);
    
    %%% FIVE-LAYER
    %induct(max(RAD)*1.2,[3.9e-3 4.1e-3],g_fl{ii},true);
    induct(max(RAD)*1.2,[2e-3 8e-3],g_fl{ii},true);
    Vfl(ii) = sum(induct(I,fils_fl{:}));
    
    % analytic solution
    thisL = defaultL(5); thisL.w = w(ii);
    Vfla(ii) = analy_sol_coil_fl(RAD,D,thisL);
end
Vfl
Vfla
save('new_pfft','w','Vhs','Vfree','Vhsa','Vfreea','Vfl','Vfla');
else
load new_pfft;
end
%%



% Get the L and R
L_c = [repmat(Vfree,size(Vhs)) real([Vhs Vfl])+Vfree];
L_a = [repmat(Vfree,size(Vhsa)) real([Vhsa Vfla])+Vfree];
R_c = -bsxfun(@times,imag([Vhs Vfl]),w(:));
R_a = -bsxfun(@times,imag([Vhsa Vfla]),w(:));

figure(1)
myplot(gca,w,L_c,w,L_a);
set(gca,'XScale','log')
xlim('auto')
ylabel('L [H]');
xlabel('w [rads^{-1}]');
legend('L_0','L_{hs}','L_{fl}');
grid on

figure(2)
myplot(gca,w,R_c,w,R_a);
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim('auto')
ylabel('\DeltaR [\Omega]');
xlabel('w [rads^{-1}]');
legend('\DeltaR_{hs}','\DeltaR_{fl}');
grid on
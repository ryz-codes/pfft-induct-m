% THIS VERSION CALCULATES AND COMPARES AGAINST ANALYTICAL SOLUTIONS
clear all
addpath([pwd '/test_fdm'])
%% Setup the frequencies and geometries to be swept
w = logspace(2,6,10);
RECALC_G = 0;

% Generate Filaments
NUM_ELE = 300; % Number of elements simulated
WID = 1e-6;
HEI = 1e-6;
D = 4e-3;
RAD = 0.1;%linspace(2.5e-2,10.5e-2,23);

% Generate a circular path
%path = genSpiralPath(RAD, -(D+HEI), NUM_ELE);
%fils = genFilsFromPath(path,WID,HEI);
%
[O L W H] = genSpiralFils( RAD, -(D+HEI), ceil(NUM_ELE/numel(RAD)), WID, HEI );
fils = {O,L,W,H}
%showFils(O,L,W,H);
I = ones(size(fils{1},1),1);

%% Setup the Green's functions
% Generate
% if RECALC_G
%     g = sweep_tph(w); % Use defaultL(1).
%     save('g_swept','g','w');
% else
    load('g_swept');
% end

%% Setup geometry
% Setup the grid and the Kernels
induct(max(RAD)*1.2,[-0.01 -0.002]);
Vfree = sum(induct(I,fils{:}));
Vfreea = analySqRing(RAD,WID);
(Vfree-Vfreea)/Vfreea

%% Iterate through each frequency
Vsub = zeros(length(w),1); % impedance calculated
Vsuba = zeros(length(w),1); % impedance analytical

for ii = 1:length(w)
    % pfft solution
    induct(max(RAD)*1.2,[-0.010 -0.002],g{ii},true);
    Vsub(ii) = sum(induct(I,fils{:}));
    
    
    % analytic solution
    thisL = defaultL(1); thisL.w = w(ii);
    Vsuba(ii) = analy_sol_coil(RAD,D,thisL);
end
Vsub
Vsuba
save('new_pfft','w','Vsub','Vfree','Vsuba');

%%

%load new_pfft;
% Get the L and R
L_c = real(Vsub);
L_a = real(Vsuba);
L_e = abs(L_c-L_a);
R_c = -imag(Vsub).*w(:);
R_a = -imag(Vsuba).*w(:);

figure(1)
myplot(gca,w,L_c,w,L_a);
set(gca,'XScale','log')
xlim('auto')
ylabel('\DeltaL [H]');
xlabel('w [rads^{-1}]');
legend('This work','Analy. sol.');
grid on

figure(2)
myplot(gca,w,R_c,w,R_a);
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim('auto')
ylabel('\DeltaR [H]');
xlabel('w [rads^{-1}]');
legend('This work','Analy. sol.');
grid on
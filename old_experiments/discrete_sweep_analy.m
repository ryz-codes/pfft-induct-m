% FREQ_SWEEP_ANALY.M
% The following script produces a plot of the analytical solution versus
% numerically simulated for a fixed discretization and for varying frequencies. 
%
% Two plots are produced, one for L0, DeltaL (half-space) and 
% DeltaL(sandwich), and one for DeltaR and DeltaR.


clear
addpath([pwd '/test_fdm'])
%% Setup the frequencies and geometries to be swept
%  Setup the Green's function if needed
n = floor(logspace(1,3,9)); % Discretization

F = [5 100]*1e3; % Frequencies to perform the analysis
RECALC_G = 1;
RECALC_Z = 1; % Re-evaluate the impedances

% Generate Filaments
NUM_ELE = 300; % Number of elements simulated
WID = 1e-6;
HEI = 1e-6;
D = 5e-3;
RAD = 0.1;%linspace(2.5e-2,10.5e-2,23);

% Generate Green's functions
if RECALC_G
    for ii = 1:length(F)
        L1 = defaultL(1);
        L1.w = 2*pi*F(ii);
        L1.bnds(3) = L1.bnds(3) *sqrt(20e3)/sqrt(F(ii));
        g_hs_d{ii} = tph(L1,0.85,0.999);
                
        L2 = defaultL(5);
        L2.w = 2*pi*F(ii);
        g_fl_d{ii} = tph(L2,0.4,0.6);     
    end
    save('g_d_swept','g_hs_d','g_fl_d');
else
    load('g_d_swept');
end

%% Calculate the ideal solutions
if RECALC_Z

% Half-space
Vhsa = zeros(1,length(F)); % impedance analytical
% Five-layer
Vfla = zeros(1,length(F)); % impedance calculated

for ii = 1:length(F)
    thisL = defaultL(1); thisL.w = 2*pi*F(ii);
    Vhsa(ii) = analy_sol_coil(RAD,D,thisL);
    
    thisL = defaultL(5); thisL.w = 2*pi*F(ii);
    Vfla(ii) = analy_sol_coil_fl(RAD,D,thisL);
end

%% Iterate through each Discretization
Vhs = zeros(length(n),length(F));
Vfl = zeros(length(n),length(F));
Vfree = zeros(length(n),1);
Vfreea = zeros(length(n),1);
for ii = 1:length(n)
    WID = 2*pi*RAD*1e-3/n(ii);
    HEI = 2*pi*RAD*1e-3/n(ii);
    
    % Freq-indep part
    % Generate circular filaments
    path = genSpiralPath(RAD, -(D+HEI/2), n(ii));
    fils_hs = genFilsFromPath(path,WID,HEI);

    path = genSpiralPath(RAD, (D+HEI/2), n(ii));
    fils_fl = genFilsFromPath(path,WID,HEI);

    % Excitation current is uniform
    I = ones(size(fils_hs{1},1),1);
    
    % Calculate the free-space component
    induct(max(RAD)*1.2,[-0.01 -0.002]);
    Vfree(ii) = sum(induct(I,fils_hs{:}));
    Vfreea(ii) = analySqRing(RAD,WID);
    
    % Freq-dependent parts
    for ij = 1:length(F)
        % Half-space
        induct(max(RAD)*1.2,[-2.5e-3 -7.5e-3],g_hs_d{ij},true);
        Vhs(ii,ij) = sum(induct(I,fils_hs{:}));

        %%% FIVE-LAYER
        induct(max(RAD)*1.2,[4.5e-3 5.5e-3],g_fl_d{ij},true);
        Vfl(ii,ij) = sum(induct(I,fils_fl{:}));
    end
end
save('disc_pfft','Vhs','Vhsa','Vfl','Vfla','Vfree','Vfreea','F','n')
else
    load('disc_pfft');
end
    
%%

% Get the L and R
L_c = real([Vhs Vfl]);
L_a = real([Vhsa Vfla]);
R_c = -bsxfun(@times,imag([Vhs Vfl]),repmat(2*pi*F,1,2));
R_a = -imag([Vhsa Vfla]).*(2*pi*[F F]);

L_0e = abs(Vfree-Vfreea)./Vfreea;
L_e = [bsxfun(@minus, L_c, L_a)];
L_e = abs(bsxfun(@rdivide, L_e, L_a));
L_e = [L_e L_0e];
R_e = bsxfun(@minus, R_c, R_a);
R_e = abs(bsxfun(@rdivide, R_e, R_a));

figure(1)
myplot(gca,n,L_e,n,L_e);
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim('auto')
ylabel('Relative Error');
xlabel('Number of Elements');
lgn_txt1 = {};
for ii = 1:length(F)
    lgn_txt1 = [lgn_txt1 {sprintf('\\DeltaL_{hs} f=%1.1fkHz',F(ii)/1000)}];
end
for ii = 1:length(F)
    lgn_txt1 = [lgn_txt1 {sprintf('\\DeltaL_{fl} f=%1.1fkHz',F(ii)/1000)}];
end
legend(lgn_txt1{:},'L_0');
grid on

figure(2)
myplot(gca,n,R_e,n,R_e);
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim('auto')
ylabel('Relative Error');
xlabel('Number of Elements');
lgn_txt2 = {};
for ii = 1:length(F)
    lgn_txt2 = [lgn_txt2 {sprintf('\\DeltaR_{hs} f=%1.1fkHz',F(ii)/1000)}];
end
for ii = 1:length(F)
    lgn_txt2 = [lgn_txt2 {sprintf('\\DeltaR_{fl} f=%1.1fkHz',F(ii)/1000)}];
end
legend(lgn_txt2);
grid on
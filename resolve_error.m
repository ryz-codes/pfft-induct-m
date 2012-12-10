

F = 5e2;
RAD= 0.1;
D=5e-3;
n = 1e3;

L1 = defaultL(1);
L1.bnds(1)=-1;
%L1.zN(1)=5000;
L1.zN(2)=2000;
L1.w = 2*pi*F;
L1.bnds(3) = L1.bnds(3) *sqrt(20e3)/sqrt(F)*5;
g_hs_d{1} = tph(L1,0.994,0.996);
%%
thisL = defaultL(1); thisL.w = 2*pi*F;
%Vhsa = analy_sol_coil(RAD,D,thisL);

%thisL = defaultL(5); thisL.w = 2*pi*F;
%Vfla = analy_sol_coil_fl(RAD,D,thisL);

WID = 2*pi*RAD*1e-3/n;
HEI = 2*pi*RAD*1e-3/n;

% Freq-indep part
% Generate circular filaments
path = genSpiralPath(RAD, -(D+HEI/2), n);
fils_hs = genFilsFromPath(path,WID,HEI);

path = genSpiralPath(RAD, (D+HEI/2), n);
fils_fl = genFilsFromPath(path,WID,HEI);

% Excitation current is uniform
I = ones(size(fils_hs{1},1),1);
%%
induct(max(RAD)*1.2,[-4.5e-3 -5.5e-3],g_hs_d{1},true);
Vhsa
Vhs = sum(induct(I,fils_hs{:}))
abs(real(Vhsa-Vhs)./real(Vhsa))
% %%% FIVE-LAYER
% induct(max(RAD)*1.2,[4.5e-3 5.5e-3],g_fl_d{1},true);
% Vfla
% Vfl = sum(induct(I,fils_fl{:}))
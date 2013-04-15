% Get calcp solution
%I = ones(N,1);
%[V] = calcp_free(I,O,L,W,H);

clear

%% Generate Filaments

NUM_ELE = 10000; % Number of elements simulated
WID = 5e-4;
HEI = 5e-4;
RAD = linspace(2.5e-2,10.5e-2,23);

% Generate a circle of filaments
[O L W H] = genCircFils( RAD, 0, ceil(NUM_ELE/numel(RAD)), WID, HEI );
I = ones(size(O,1),1);
%showFils(O,L,W,H);
sum(resist(I,O,L,W,H,0.5));

%% Setup pfft
induct(O,O+L+W+H);
%%
tic
Vpfft = induct(I,O,L,W,H);
toc

fprintf('\n\n');
    fprintf('P-FFT result: %g\n',sum(Vpfft));

%% Run direct solve
if NUM_ELE <= 2e3
    tic
    Vcalcp = calcp_free(I,O,L,W,H);
    toc
    
    fprintf('Calcp result: %g\n',sum(Vcalcp));
    fprintf('Peak error: %g\n',max(abs(Vpfft-Vcalcp)./Vcalcp));
    fprintf('Mean error: %g\n',mean(abs(Vpfft-Vcalcp)./Vcalcp));
end

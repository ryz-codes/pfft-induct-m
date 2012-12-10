% This version does not set up the inv_r portion

%clear
%% Generate Filaments

NUM_ELE = 100; % Number of elements simulated
WID = 4e-4;
HEI = 4e-4;
D = 3e-3;
RAD = 0.1%linspace(2.5e-2,10.5e-2,23);

% Generate a circle of filaments
[O L W H] = genSpiralFils( RAD, -(D+HEI), ceil(NUM_ELE/numel(RAD)), WID, HEI );
showFils(O,L,W,H);
I = ones(size(O,1),1);

%% Setup pfft
load default1.mat;
induct(max(RAD)*1.2,[-0.010 -0.002],g,true);

%%
Vpfft = induct(I,O,L,W,H);
fprintf('\n\n');
fprintf('P-FFT result: %g + j%g\n',real(sum(Vpfft)),imag(sum(Vpfft)));
    
%% Run direct solve
if NUM_ELE < 5e2
    tic
    Vcalcp = calcp(I,O,L,W,H);
    toc
end

%%
if ~isempty(Vcalcp)
    fprintf('Calcp result: %g + j%g\n',real(sum(Vcalcp)),imag(sum(Vcalcp)));
    fprintf('Peak error: %g\n',max(abs(Vpfft-Vcalcp)./abs(Vcalcp)));
    fprintf('Mean error: %g\n',mean(abs(Vpfft-Vcalcp)./abs(Vcalcp)));
end
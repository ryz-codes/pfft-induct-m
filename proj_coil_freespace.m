% Get calcp solution
%I = ones(N,1);
%[V] = calcp_free(I,O,L,W,H);

clear 

%% Generate Filaments
REFINE = 4;
NUM_ELE = 4000;
fils = genCoilFils(REFINE,NUM_ELE);

% Generate a circle of filaments
I = ones(size(fils{1},1),1)/REFINE^2;
sum(resist(I,fils{:},1));

%% Setup pfft
induct(0.123,[-0.028 0.006]);
%%
tic
Vpfft = induct(I,fils{:});
toc

fprintf('\n\n');
    fprintf('P-FFT result: %g\n',sum(Vpfft));

%% Run direct solve
if NUM_ELE <= 600
    tic
    Vcalcp = calcp_free(I,fils{:});
    toc
    
    fprintf('Calcp result: %g\n',sum(Vcalcp));
    fprintf('Peak error: %g\n',max(abs(Vpfft-Vcalcp)./Vcalcp));
    fprintf('Mean error: %g\n',mean(abs(Vpfft-Vcalcp)./Vcalcp));
end

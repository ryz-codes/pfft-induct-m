% FREQ_SWEEP_PROJ
GREEN_FN = 0;
FREE_SPACE = 0;
D_OFF = (0.709-0.295/2)*25.4e-3; % dist programmed into our green's functions
NEW_D = (0.106+0.295/2)*25.4e-3%(0.11+0.295*0.5)*25.4e-3;
% F = [0.1; 1; 10; 100]*1e3;
F= [0.1
    0.2
    0.5
    1
    2
    5
    10
    20
    50
    100
    200
    500
    1000]*1e3;

%% Generate Filaments
REFINE = 4;
NUM_ELE = 4000;
fils = genCoilFils(REFINE,NUM_ELE,D_OFF-NEW_D);

I = ones(size(fils{1},1),1)/REFINE^2;
sum(resist(I,fils{:},0.6));

%% Free-Space case
if FREE_SPACE
induct(0.123,[-0.028 0.006+D_OFF-NEW_D]);
tic
Vpfft = induct(I,fils{:});
toc

fprintf('P-FFT result: %g\n',sum(Vpfft));
end

%% One plate case
if GREEN_FN
    clear g
    for ii = 1:length(F)
        L = defaultL(10); % Plate case
        L.w = 2*pi*F(ii);
        L
        g{ii} = tph(L,0.85,0.999);
    end
     save('proj_oneplate','F','g');
else
    load('proj_oneplate');
end
%%
clear Vsub
for ii = 1:length(F)
    induct(0.126,[-0.028 D_OFF-4e-3],g{ii},true);
    Vsub(ii) = sum(induct(I,fils{:}))
end
%%
disp(-real(Vsub).'*1e6)
disp(-imag(Vsub).'.*(2*pi*F))
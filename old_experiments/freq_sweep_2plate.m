% This is the two-plate case
%
% DIMENSIONS used
% The top of the bottom plate to the bottom of the coil:
%   1.0825 inch +- 0.01
% The bottom of the coil to the bottom of the top plate without spacers:
%   0.094 inch +- 0.015
% Thickness of spacers
%   0.7095 inch +- 0.0005
% Top of the coil to the bottom of the top plate without spacers:
%   0.11 inch
% Height of the coil
%   0.295 inch

% FREQ_SWEEP_PROJ
GREEN_FN = 0;
FREE_SPACE = 0;
D_OFF = 0; % dist programmed into our green's functions
NEW_D = (0.709-0.295)*25.4e-3;%(0.106+0.295/2)*25.4e-3%(0.11+0.295*0.5)*25.4e-3;
 F = [0.1; 1; 10; 100]*1e3;
%F= [0.1
%     0.2
%     0.5
%     1
%     2
%     5
%     10
%     20
%     50
%     100
%     200
%     500
%     1000]*1e3;

%% Generate Filaments
REFINE = 4;
NUM_ELE = 4000;
fils = genCoilFils(REFINE,NUM_ELE,-NEW_D-(0.295/2*25.4e-3));
uLz = fils{2};
uLz = uLz(:,3)./sqrt(sum(fils{2}.^2,2));
ut = sqrt(1- uLz.^2); % transverse component

I = ut/REFINE^2;
%I = ones(size(fils{1},1),1)/REFINE^2;
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
        L = defaultL(11); % 2-Plate case
        L.w = 2*pi*F(ii);
        L
        g{ii} = tph(L,0.05,0.95);
    end
     save('proj_twoplate','F','g');
else
    load('proj_twoplate');
end
%%
clear Vsub

% figure out the bounding box.
% upper bound
zup = max(fils{1}+fils{4});
zup = zup(3);
zdown = min(fils{1});
zdown = zdown(3);
dz = (zup-zdown)/(2^8-2);
for ii = 1:length(F)
    induct(0.126,[zdown-dz zup+dz],g{ii},true);
    a = induct(I,fils{:});
    Vsub(ii) = sum(a)
end
%%
disp(-real(Vsub).'*1e6)
disp(-imag(Vsub).'.*(2*pi*F))
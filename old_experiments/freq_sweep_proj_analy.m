% Generates the results using analytical formulas
%--------------------------------------------------------------------------
% PHYSICAL MEASUREMENTS
% Lisa's measurements
% first set of turns
a1(1)= 1+15/16-14/32; %inch
a1(2)= 1+15/16;
a1 = a1*25.4;
N1 = 7;

% second set of turns
a2(1) = 3.25-18/32;
a2(2) = 3.25;
a2 = a2*25.4;
N2 = 9;

% third set of turns
a3(1) = 4.5-22/32;
a3(2) = 4.5;
a3 = a3*25.4;
N3 = 12;
 
a = [linspace(a1(1),a1(2),N1) linspace(a2(1),a2(2),N2)...
    linspace(a3(1),a3(2),N3)]*1e-3;

% the wires themselves
wid = (0.063*25.4)*1e-3; %[m] 0.063 in +- 0.002in
hei = 7.493e-3; %[m] 0.295 in +- 0.005in

dist_board = 0.709*25.4e-3; % dist board to board
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
z = dist_board - hei/2;
z = 0.106*25.4e-3+ hei/2;
%%
L0 = analy_coils_free(a,mean([wid, hei]))
%% Calculate
L = defaultL(10);
for ii = 1:length(F)
    L.w = 2*pi*F(ii);
    V(ii) = analy_coils_plate(a,z,L);
    ii
end

%% Process
DL = -real(V).'*1e6
DR = -imag(V).'.*F(:)*2*pi
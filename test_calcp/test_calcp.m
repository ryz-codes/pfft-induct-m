clc

% Test circle of filaments
R = 0.2;
Z = -0.02;
w = 20e3 *2 *pi;

N = 100;
wid = 15e-3;
hei = 15e-3;


[O L W H] = genCircFils( R, Z, N, wid, hei );
showFils(O,L,W,H);

%%

I = ones(N,1);
[V M] = calcp(I,O,L,W,H);
disp('Calculated');
M1 = sum(V)
disp('Ideal solution');
M2 = analy_sol_hs(R,Z,Z,w)
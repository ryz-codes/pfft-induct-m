scale = 1e-2;

O = [0 1 0 0 1.2 0]';
L = [1 0 0 1 0 0]'*scale;
W = [0 1 0 0 1 0]'*scale;
H = [0 0 1 0 0 1]'*scale;

N=2e6;
O = kron(ones(1,N),O);
L = kron(ones(1,N),L);
W = kron(ones(1,N),W);
H = kron(ones(1,N),H);
tic
mut(O,L,W,H);
toc


%%

NUM_ELE = 10; % Number of elements simulated
WID = 0.01;
HEI = 0.01;
RAD = 1;
[O L W H] = genCircFils( RAD, 0, NUM_ELE, WID, HEI );
I = ones(NUM_ELE,1);
tic
Vcalcp = calcp_free(I,O,L,W,H);
toc
b = sum(Vcalcp);
    fprintf('Calcp solution: %g\n',b);
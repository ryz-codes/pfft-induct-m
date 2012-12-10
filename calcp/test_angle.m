% To test the accuracy of Matt's brick to brick mutual inductance code for
% cases where the bricks are at an angle to each other. My suspicion is
% that the quadrature rule used is not very good for close-by interactions.

WIDTH = 0.1;
HEIGHT = 0.1;
LENGTH = 1;
N = 1;

kp = @(a) [kron(a,ones(size(a,1),1)) kron(ones(size(a,1),1),a)];

% First set
vert = linspace(0,1,N+1).';
O1 = -vert(1:N);
L1 = -vert(2:N+1)-O1;
O1 = [O1 zeros(N,2)];
L1 = [L1 zeros(N,2)];
W1 = repmat([0 1 0]*WIDTH,N,1);
H1 = repmat([0 0 1]*HEIGHT,N,1);

O1 = O1 - W1/2 - H1/2;

vert2 = bsxfun(@times,vert,1/sqrt(2)*[1 -1 0]);
O2 = vert2(1:N,:);
L2 = vert2(2:N+1,:)-O2;
W2 = repmat(1/sqrt(2)*[-1 -1 0]*WIDTH,N,1);
H2 = repmat([0 0 1]*HEIGHT,N,1);

O2 = O2 - W2/2 - H2/2;

O = [O1;O2]; L = [L1;L2]; W = [W1;W2]; H = [H1;H2];
showFils(O,L,W,H)

V = mutual(kp(O)',kp(L)',kp(W)',kp(H)');
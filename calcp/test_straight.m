WIDTH = 0.1;
HEIGHT = 0.1;
LENGTH = 1;
N = 10;

kp = @(a) [kron(a,ones(N,1)) kron(ones(N,1),a)];

vert = linspace(0,1,N+1);
O = vert(1:N).';
L = vert(2:N+1).'-O;
O = [O zeros(N,2)];
L = [L zeros(N,2)];
W = repmat([0 1 0]*WIDTH,N,1);
H = repmat([0 0 1]*HEIGHT,N,1);

O = O - W/2 - H/2;
%showFils(O,L,W,H)

V = mutual(kp(O)',kp(L)',kp(W)',kp(H)');
disp(sum(V))
%%
V2 = filquad(kp(O),kp(L),kp(W),kp(H));

%%
% Extract mutual terms
ind = eye(N);
ind = ~logical(ind(:));

V= V(ind);
V2 = V2(ind);

%%
disp('V | V2 | V-V2 | rel(V-V2)');
disp([V V2 V-V2])
disp((V-V2)./V)
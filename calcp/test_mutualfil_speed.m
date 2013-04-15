% a = [0     0     0];
% b = [0.0100         0         0];
% c = [0    0.0100         0];
% d = [0.0100    0.0100         0];
% 
% N = 1e6;
% 
% a = rand(N,6);
% 
% % tic
% % for ii= 1:N
% %     mutualfil(a(ii,1),a(ii,2),a(ii,3),a(ii,4),a(ii,5),a(ii,6));
% % end
% % toc
% tic
% mutualfil(a(:,1),a(:,2),a(:,3),a(:,4),a(:,5),a(:,6));
% toc

%%

scale = 3e-3;

O = [0 0 0 0 0 0]';
L = [1 0 0 1 0 0]';
W = [0 1 0 0 1 0]'*scale;
H = [0 0 1 0 0 1]'*scale;
%showFils(O',L',W',H')
N=1;
O = kron(ones(1,N),O);
L = kron(ones(1,N),L);
W = kron(ones(1,N),W);
H = kron(ones(1,N),H);
tic
V1 = mutual(O,L,W,H);
toc
V1(1)
%tic
%V2 = mutual2(O,L,W,H);
%toc
%disp(mean(V1-V2));


%%
% Test with random inputs
% N = 100000
% O = randi(100, 6,N);
% L = randi(100, 6,N);
% W = randi(100, 6,N);
% H = randi(100, 6,N);
% tic
% mutual2(O,L,W,H);
% toc
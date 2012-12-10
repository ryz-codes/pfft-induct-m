function [Hr Tr r zh zt] = tableTPH(L,lbnd,ubnd)
%TABLETPH Given the layer structure L, makes a T+H lookup table for the
%layer selected by the flag coil_layer.
% SYSTEM CONSTANTS
INTERP_COLS = 10;
if nargin == 1
    lbnd = 0.01;
    ubnd = 0.99;
end

% Identify coil layer
coil_layer = L.coil_layer;     
zN = L.zN(coil_layer);

% Assert that this layer can contain the coil
assert(L.sig(coil_layer)==0 &&...
    L.mu_r(coil_layer)==1, 'Coil layer must be free-space!');

% Identify upper and lower boundaries
bnd_a = L.bnds(coil_layer+1);
bnd_b = L.bnds(coil_layer);

% Make the "box" 95% of the way to either boundary.
z_val = linspace(bnd_b,bnd_a,zN);
z = z_val(round(lbnd*zN):round(ubnd*zN));

% Distribute runs evenly between two ends
N = length(z);
ip = round(linspace(1,N,INTERP_COLS));

tstart = tic;
A = cell(INTERP_COLS,1);
for ii = 1:INTERP_COLS
    fprintf('Run %u/%u: z=%g... ',ii,INTERP_COLS,z(ip(ii))); tic;
    [A{ii},r,zout] = edd2([z(1),z(N)],z(ip(ii)),L);
    fprintf('Ok (%gs)\n',toc);
end
tFDM = toc(tstart);


% Do each r value
Tr = [];
Hr = [];
tic;
for ii = 1:length(r)
    fprintf('Split %u/%u: r=%g... ',ii,length(r),r(ii)); 
     G = zeros(N,INTERP_COLS);
     for ij = 1:INTERP_COLS
         G(:,ij) = (A{ij}(ii,:)).';
     end
    
    [Tr1 Hr1 FLAG RELRES ITER] = splitTPH(G,ip);
    Tr = [Tr Tr1];
    Hr = [Hr Hr1];
    fprintf('%u iters (res=%g)\n',ITER,RELRES);
end

disp('====================');
fprintf('Total FDM time=%gs\n',tFDM);
fprintf('Total split time=%gs\n',toc);

% Prepare the guide vectors.
r(1) = 0; % smudge the first r point to make it zero
zh = linspace(2*z(1),2*z(N),2*N-1);
zt = linspace(0,z(N)-z(1),N);

end


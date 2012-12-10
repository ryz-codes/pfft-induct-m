zp = -5e-3;


%%
rN = 60;
rbnd = 1; % r upper bound.

zN1 = 1800;
zN2 = 500;

w = 2*pi*20e3;
mu = 4*pi*1e-7;
sig = 8e6;
psn = -1j*w*mu*130*sig;
skin_d = sqrt(2/w/sig/mu/130);

dz1 = 0.5/1800;
mu_r1 = 1;
dz2 = skin_d/100;
mu_r2 = 130;

tic
%%
% Begin trying to piece together the boundary condition
[A1 r] = blk(rN,rbnd,zN1,dz1);
A2 = blk(rN,rbnd,zN2,dz2,psn);
toc
tic
rows = zN1+[0,1,1,2,2,2,2,3];
cols = zN1+[1,1,2,0,1,2,3,2];
ents = [1/dz1^2, ...
        -1,1, ...
        [1,-1]./mu_r1./dz1 , [-1,1]./mu_r2./dz2, ...
        1/dz2^2];
Abnd = sparse(rows,cols,ents,zN1+zN2+2,zN1+zN2+2);
Abnd = kron(Abnd, speye(rN));

A = blkdiag(A1, sparse(rN*2,rN*2), A2);
A = A + Abnd;

%%
b1 = mu./4./pi./sqrt(r.^2 + (zp).^2);
b2 = mu./4./pi.*(zp)./(sqrt(r.^2 + (zp).^2).^3);

rows = rN*zN1+(1:(2*rN));

b = sparse(rows,1,[b1;b2],length(A),1);

%%
ind1 = zN1+1;
x = A \ b;
x = reshape(x,rN,[]);
xout = x(:,ind1);

%%
[xr r2] = edd_fast(0,-5e-3);
%%
figure;plot(r,real(xout),r2,real(xr));
figure;plot(r,imag(xout),r2,imag(xr));
xlim([0 0.1]);
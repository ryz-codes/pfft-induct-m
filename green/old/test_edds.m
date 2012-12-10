
zp = -49e-3;
z = -49e-3;
%%
[x1 r1 z1] = edd_fast(z,zp);

%%
L = defaultL(4);
tic
[x2 r2 z2] = edd2(z,zp,L);
toc

%%
figure(1);
subplot(211);
plot(r1,real(x1),r2,real(x2));
xlim([0 0.1]);
subplot(212);
plot(r1,imag(x1),r2,imag(x2));
xlim([0 0.1]);
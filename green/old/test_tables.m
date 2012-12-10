%% Setup tables
clear;
L = defaultL(4);
[H T ri zh zt] = tableTPH(L);

%% Lookup
d=-4e-3;
z=-4e-3;
r=0:0.001:0.1;

[Gr] = lookupTPH(ri,zh,zt,H,T,r,z,d);
[Gfdm r2] = edd_fast(z,d);
%[Gfdm r2] = edd2(z,d,L);


figure (2);
clf(2);
subplot(211);
hold on
plot(r,real(Gr),'o')
plot(r2,real(Gfdm))
xlim([0,0.1])

subplot(212);
hold on
plot(r,imag(Gr),'o')
plot(r2,imag(Gfdm))
xlim([0,0.1])

%% Lookup
d=-40e-3;
z=linspace(-100e-3,-1e-3,100);
r=0;

[Gr] = lookupTPH(ri,zh,zt,H,T,r,z,d);
[Gfdm r2 zout] = edd2([z(1),z(end)],d,L);


figure (2);
clf(2);
subplot(211);
hold on
plot(z,real(Gr),'o')
plot(zout,real(Gfdm(1,:)))
%xlim([0,0.1])

subplot(212);
hold on
plot(z,imag(Gr),'o')
plot(zout,imag(Gfdm(1,:)))
%xlim([0,0.1])



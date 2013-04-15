%% Setup tables
clear;
load('default4.mat');

%% Lookup
d=-3e-3;
z=-3e-3;
r=0:0.001:0.1;
w = 20e3*2*pi;

[Gr] = g.lookup(r,z,d);
[Gfdm r2] = analy_sol_hs(z+d,w);

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


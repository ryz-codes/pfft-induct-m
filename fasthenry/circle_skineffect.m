% Specify the geometry
N = 100; % Number of segments
rad = 2.5e-4; % Cross-sectional radius
R = linspace(0.025,0.105,2); % Overall radius of circle
f = logspace(1,6,20);

% Make the path
verts = genSpiralPath( R, 0, N);

% Make the filaments
geom = struct;
[geom.O, geom.L, geom.W, geom.H] = ...
    path2segs(verts,rad);

% Run FastHenry
[Zt] = fasthenry(geom,f);

% Plot cross-sectional current density
figure(2);
subplot(121);semilogx(f,imag(Zt)./(2*pi*f)); title('Inductance vs f');
subplot(122);loglog(f,real(Zt)); title('Resistance vs f');
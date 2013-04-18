% Specify the geometry
N = 1000; % Number of segments
numlitz = 15;
bundle_rad = 7.493e-3/2/4; %
f = [0.1,1,3,5,10,20,30,50,100,200,300,500,1000]*1e3;

% the wires themselves
wid = (0.063*25.4)*1e-3; %[m] 0.063 in +- 0.002in
hei = 7.493e-3; %[m] 0.295 in +- 0.005in
filfac = 0.6;
rad = sqrt(wid*hei*filfac/pi/numlitz); % Do cross-sectional radius to maintain DC resistance.

% Make the path
fils = genCoilFils(1,N);
verts = fils{1};
%verts = genSpiralPath( [0.025,0.105], 0, 50);
%%
litz_verts = path2litz(verts,bundle_rad,numlitz,numlitz);
length(verts)

% Make the filaments
geom = struct;
geom.O = cell(numlitz,1);
geom.L = cell(numlitz,1);
geom.W = cell(numlitz,1);
geom.H = cell(numlitz,1);
for ii = 1:numlitz
[t1,t2,t3,t4] = path2segs(litz_verts{ii},rad);
    geom.O{ii} = t1;
    geom.L{ii} = t2;
    geom.W{ii} = t3;
    geom.H{ii} = t4;
end
%%
%showFils(vertcat(geom.O{1}{:}),vertcat(geom.L{1}{:}),vertcat(geom.W{1}{:}),vertcat(geom.H{1}{:}));
%axis equal
%geom2.O = geom.O{1}; geom2.L = geom.L{1}; geom2.W = geom.W{1}; geom2.H = geom.H{1};
%%
% Run FastHenry
t=tic
[Zt] = fasthenry(geom,f);
toc(t)
%% Plot cross-sectional current density
figure(1);
Rt = real(Zt);
Lt = imag(Zt)./(2*pi*f);
subplot(121);semilogx(f,imag(Zt)./(2*pi*f)); title('Inductance vs f');
subplot(122);loglog(f,real(Zt)); title('Resistance vs f');
save Zt f Zt Rt Lt;
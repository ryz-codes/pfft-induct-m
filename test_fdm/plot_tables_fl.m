% The following file plots the subtracted Green's function for a variety of
% z', with each normalized to Gsub(0,0,z');
% A number of hacks were used to make this following code work for
% defaultL(5).

clear
%% Setup tables
RMAX = 0.1;
zp= [0.01 0.1 3 5]*1e-3;
w = 10e3*2*pi;

% Make text for legend
lgn_txt = cell(size(zp));
for ii = 1:length(zp)
    lgn_txt{ii} = sprintf('z''=%gmm',abs(zp(ii))*1e3);
end

% Lookup
L = defaultL(5);
Gfdm_r = cell(size(zp));
Gfdm_z = cell(size(zp));

zlims = [L.bnds(L.coil_layer) L.bnds(L.coil_layer+1)];

L.w = w;
for ii = 1:length(zp)
    % Run FDM
    fprintf('Running FDM with z''= %g... ',zp(ii));
    [fdm_res r z] = fdm_run(zlims,zp(ii),L);
    fprintf('DONE\n');
    
    % Extract result for G vs r plot
    z_idx = find(z<=0,1,'last'); % This is the first index where z=0 in the coil domain
    Gfdm_r{ii} = fdm_res(:,z_idx);
    
    % Extract result for G vs z plot
    Gfdm_z{ii} = fdm_res(1,:).';
    
    % Analytical solution vs r
    [Ganaly{ii} r2] = analy_sol_fl(0,zp(ii),L);
    Ganaly{ii} = Ganaly{ii}.';
end


%%
% Transform cell into matrix
if iscell(Gfdm_r)
    Gfdm_r = cell2mat(Gfdm_r);
    Gfdm_z = cell2mat(Gfdm_z);
    Ganaly = cell2mat(Ganaly);
end

% Plot G vs r
figure(1);
clf(1);
h = gca;
[G2 G1] = my_norm(abs(Ganaly),abs(Gfdm_r));
myplot(h,r,G1,r2,G2)
xlim(h,[0,0.05])
ylim(h,[0, 1]);
grid on;
legend(lgn_txt);
title('Normalized ||G_{sub}(\rho,z=0,z'')||');
xlabel('\rho [m]');
%%
% Plot G vs z
figure(2);
clf(2);
h = gca;
[G2 G1] = my_norm(abs(Gfdm_z),abs(Gfdm_z));
z1 = downsample(z,20);
G1a = downsample(G1,20);
myplot(h,z1,G1a,z,G2)
xlim(h,[0,0.01])
ylim(h,[0, 1.5]);
grid on;
legend(lgn_txt);
title('Normalized ||G_{sub}(\rho=0,z,z'')||');
xlabel('z [m]');
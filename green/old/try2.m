clear
zp1 = -5e-2; % Lowest point
zp2 = -5e-3; % Highest point

% Collect data from the first two runs
[Abot,r,z] = edd_fast(zp1,[zp1,zp2]);
[Atop] = edd_fast(zp2,[zp1,zp2]);

% Figure out where to get other two runs
N = length(z);
ip2 = round(N/3);
ip3 = round(N*2/3);

A2 = edd_fast(z(ip2),[zp1,zp2]);
A3 = edd_fast(z(ip3),[zp1,zp2]);

% Combine
ip = [1,ip2,ip3,N];

% Do each r value
Tr = [];
Hr = [];
for ii = 1:length(r)
    G = [Abot(ii,:);A2(ii,:);A3(ii,:);Atop(ii,:)].';
    [Tr1 Hr1] = splitTPH(G,ip);
    Tr = [Tr Tr1];
    Hr = [Hr Hr1];
end
%%
zh = linspace(2*z(1),2*z(N),2*N-1);
zt = linspace(0,z(N)-z(1),N);


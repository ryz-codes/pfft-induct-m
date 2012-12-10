function [Hr Tr r zh zt] = init_table(zp1,zp2)
if nargin < 2 % reset to defaults
    zp1 = -5e-2; % Lowest point
    zp2 = -5e-3; % Highest point
end

% Collect data from the first two runs
[Abot,r,z] = edd_fast([zp1,zp2],zp1);
[Atop] = edd_fast([zp1,zp2],z(end));

% Figure out where to get other four runs
N = length(z);

ip2 = round(N/5);
ip3 = round(N*2/5);
ip4 = round(N*3/5);
ip5 = round(N*4/5);

A2 = edd_fast([zp1,zp2],z(ip2));
A3 = edd_fast([zp1,zp2],z(ip3));
A4 = edd_fast([zp1,zp2],z(ip4));
A5 = edd_fast([zp1,zp2],z(ip5));

% Combine
ip = [1,ip2,ip3,ip4,ip5,N];

% Do each r value
Tr = [];
Hr = [];
for ii = 1:length(r)
    G = [Abot(ii,:);A2(ii,:);A3(ii,:);A4(ii,:);A5(ii,:);Atop(ii,:)].';
    [Tr1 Hr1] = splitTPH(G,ip);
    Tr = [Tr Tr1];
    Hr = [Hr Hr1];
end

% smudge the first r point to make it zero
r(1) = 0;

%%
zh = linspace(2*z(1),2*z(N),2*N-1);
zt = linspace(0,z(N)-z(1),N);


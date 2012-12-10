function [O L W H] = genSpiralFils( R, Z, N, wid, hei )
%genCircFils Generates a spiral of filaments
%   Detailed explanation goes here

num_turns = numel(R);
if num_turns == 1 % Make a circle
    [O L W H] = genCircFils(R,Z,N,wid,hei);
    return
end

N_tot = N*num_turns;
r_path = linspace(min(R),max(R),N_tot+1)';
th_path = linspace(0,2*pi*num_turns,N_tot+1)';
verts = [r_path.*cos(th_path), r_path.*sin(th_path), repmat(Z,N_tot+1,1)];
clear r_path th_path

O = verts(1:N_tot,:);
L = verts(2:N_tot+1,:) - O;
ul = bsxfun(@rdivide,L, (sqrt(sum(L.^2,2))));
W = ul * [0 -1 0; 1 0 0; 0 0 1] * wid; % rotate
H = repmat([0 0 1],N_tot,1) * hei; % copy

% Shift everything down W and H
O = O - 0.5*W - 0.5*H;
    
end

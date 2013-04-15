function [verts] = genSpiralPath( R, Z, N)
%genCircFils Generates a spiral of filaments
%   Detailed explanation goes here

num_turns = numel(R);

N_tot = N*num_turns;
r_path = linspace(min(R),max(R),N_tot+1)';
th_path = linspace(0,2*pi*num_turns,N_tot+1)';
verts = [r_path.*cos(th_path), r_path.*sin(th_path), repmat(Z,N_tot+1,1)];
    
end

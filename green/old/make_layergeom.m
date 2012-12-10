function [L] = make_layergeom()
%MAKE_LAYERGEOM Summary of this function goes here
%   Detailed explanation goes here
    L = struct;
    
    % THE FOLLOWING DEFAULT VALUES IMPLEMENT THE STANDARD HALF-SPACE
    % PROBLEM WITH THE COIL LOCATED Z<0 AND THE INTERFACT LOCATED AT Z=0.
    % MU_R AND SIGMA USED HERE CORRESPOND TO THOSE USED BY ACERO.
    
    % Number of layers
    L.num = 2;
    
    % Coil layer is at this layer. Layers are numbered from low z to high
    % z.
    L.coil = 1;

    % Relative mu of each layer, starting at layer with lowest z
    L.mu_r = [1 130];
    
    % Conductivity of each layer
    L.sig = [0 8e6];
    
    % Boundary locations in z
    L.bound = [0];
    
    % Excitation using "above coil" equations
    L.excite_above = 1;
    
    % Excitation using "below coil" equations
    L.excite_below = 0;
    
    
    

end


function [L] = defaultL(SELECT)
%DEFAULTL Summary of this function goes here
%   Detailed explanation goes here
% The layer geometry is defined such that the bottom of the first layer
% starts at z=offset. 
L = struct;

% General parameters
L.w = 2*pi*2e4;
L.mu = 4*pi*1e-7;

if SELECT ==1;
    %---------------------------------------------------------------------- 
    %  Coil below interface
    %----------------------------------------------------------------------  

    % Define region
    L.layerN = 2;

    % Layer 1 (approaching from -inf)
    L.names{1} = 'Air containing coil';
    L.mu_r(1) = 1;
    L.zN(1) = 1800;
    L.bnds(1) = -0.5; % boundary below
    L.sig(1) = 0;

    % Layer 2
    L.names{2} = 'Cast Iron (from Acero)';
    L.mu_r(2) = 130;
    L.zN(2) = 500;
    L.sig(2) = 8e6;
    skin_d = sqrt(2/L.w/L.sig(2)/L.mu_r(2)/L.mu);
    L.bnds(2) = 0; 

    % Upper bound
    L.bnds(3) = 5*skin_d;

    % Excitation
    L.coil_layer = 1;
end
if SELECT == 2
    %---------------------------------------------------------------------- 
    %  Coil above interface
    %----------------------------------------------------------------------

    % Define region
    L.layerN = 2;

    % Layer 1 (approaching from -inf)
    L.names{1} = 'Cast Iron (from Acero)';
    L.mu_r(1) = 130;
    L.zN(1) = 200;
    L.sig(1) = 8e6;
    skin_d = sqrt(2/L.w/L.sig(1)/L.mu_r(1)/L.mu);
    L.bnds(1) = -10*skin_d;
    
    % Layer 2
    L.names{2} = 'Air containing coil';
    L.mu_r(2) = 1;
    L.zN(2) = 200;
    L.sig(2) = 0;    
    L.bnds(2) = 0; 

    % Upper bound
    L.bnds(3) = 10;

    % Excitation
    L.coil_layer = 2;
end

if SELECT == 3
    %---------------------------------------------------------------------- 
    %  Coil above interface with coarse mesh above air layer
    %----------------------------------------------------------------------

    % Define region
    L.layerN = 3;

    % Layer 1 (approaching from -inf)
    L.names{1} = 'Cast Iron (from Acero)';
    L.mu_r(1) = 130;
    L.zN(1) = 500;
    L.sig(1) = 8e6;
    skin_d = sqrt(2/L.w/L.sig(1)/L.mu_r(1)/L.mu);
    L.bnds(1) = -5*skin_d;
    
    % Layer 2
    L.names{2} = 'Air containing coil';
    L.mu_r(2) = 1;
    L.zN(2) = 1000;
    L.sig(2) = 0;    
    L.bnds(2) = 0; 
    
    % Layer 3
    L.names{3} = 'Coarsely meshed air';
    L.mu_r(3) = 1;
    L.zN(3) = 500;
    L.sig(3) = 0;    
    L.bnds(3) = 0.05; 

    % Upper bound
    L.bnds(4) = 0.5;

    % Excitation
    L.coil_layer = 2;
end

if SELECT ==4;
    %---------------------------------------------------------------------- 
    %  Coil below interface
    %----------------------------------------------------------------------  

    % Define region
    L.layerN = 4;

    % Layer 1 (approaching from -inf)
    L.names{1} = 'V Coarsely meshed air';
    L.mu_r(1) = 1;
    L.zN(1) = 100;
    L.bnds(1) = -1; % boundary below
    L.sig(1) = 0;
    
    % Layer 1 (approaching from -inf)
    L.names{2} = 'Coarsely meshed air';
    L.mu_r(2) = 1;
    L.zN(2) = 200;
    L.bnds(2) = -0.3; % boundary below
    L.sig(2) = 0;
    
    % Layer 2
    L.names{3} = 'Air containing coil';
    L.mu_r(3) = 1;
    L.zN(3) = 300;
    L.bnds(3) = -0.1; % boundary below
    L.sig(3) = 0;

    % Layer 3
    L.names{4} = 'Cast Iron (from Acero)';
    L.mu_r(4) = 130;
    L.zN(4) = 300;
    L.sig(4) = 8e6;
    skin_d = sqrt(2/L.w/L.sig(4)/L.mu_r(4)/L.mu);
    L.bnds(4) = 0; 

    % Upper bound
    L.bnds(5) = 5*skin_d;

    % Excitation
    L.coil_layer = 3;
end

if SELECT ==5;
    %---------------------------------------------------------------------- 
    %  Five layer, sandwich configuration
    %----------------------------------------------------------------------  
    PLATE_THICK = 1e-3; % Skin depth of aluminum = 0.3mm at 1MHz
    SEP_DIST = 0.01;
    
    % Define region
    L.layerN = 5;

    % Layer 1 (approaching from -inf)
    L.names{1} = 'Air';
    L.mu_r(1) = 1;
    L.zN(1) = 100;
    L.bnds(1) = -10; % boundary below
    L.sig(1) = 0;
    
    % Layer 2 (approaching from -inf)
    L.names{2} = 'Aluminium';
    L.mu_r(2) = 1;
    L.zN(2) = 100;
    L.sig(2) = 3.5e7;
    L.bnds(2) = -(PLATE_THICK); % boundary below
    
    % Layer 3
    L.names{3} = 'Air containing coil';
    L.mu_r(3) = 1;
    L.zN(3) = 200;
    L.bnds(3) = 0; % boundary below
    L.sig(3) = 0;

    % Layer 4
    L.names{4} = 'Aluminium';
    L.mu_r(4) = 1;
    L.zN(4) = 100;
    L.sig(4) = 3.5e7;
    L.bnds(4) = SEP_DIST; 
    
    % Layer 5
    L.names{5} = 'Air';
    L.mu_r(5) = 1;
    L.zN(5) = 100;
    L.bnds(5) = SEP_DIST+PLATE_THICK; % boundary below
    L.sig(5) = 0;
    
    % Upper bound
    L.bnds(6) = 10;

    % Excitation
    L.coil_layer = 3;
end
if SELECT ==10;
    % PROJECT CONFIGURATION
    %---------------------------------------------------------------------- 
    %  Coil below interface
    %----------------------------------------------------------------------  
    PLATE_THICK = 25.4/8*1e-3; % 1/8 inch thick (0.125 +- 0.0005 in, 0.4%)
    % Spacers 0.709 in +- 0.0005, coil thickness 0.29-0.3 in
    % The following distance is between the coil mid-point and the plates
    COIL_PLATE_SEP_DIST = (0.7095-0.295/2)*25.4e-3;
    L.w = 2*pi*1e3;
    
    % Define region
    L.layerN = 3;

    % Layer 1 (approaching from -inf)
    L.names{1} = 'Air containing coil';
    L.mu_r(1) = 1;
    L.zN(1) = 1800;
    L.bnds(1) = -0.5; % boundary below
    L.sig(1) = 0;

    % Layer 2
    L.names{2} = 'Annealed Copper';
    L.mu_r(2) = 1;
    L.zN(2) = 500;
    L.sig(2) = 5.8e7;
    L.bnds(2) = COIL_PLATE_SEP_DIST; 

    % Upper bound
    L.names{3} = 'Air';
    L.mu_r(3) = 1;
    L.zN(3) = 200;
    L.sig(3) = 0;
    L.bnds(3) = COIL_PLATE_SEP_DIST+PLATE_THICK;
    
    L.bnds(4) = 0.5;

    % Excitation
    L.coil_layer = 1;
end

if SELECT ==11;
    % PROJECT CONFIGURATION
    %---------------------------------------------------------------------- 
    %  Five layer, sandwich configuration
    %----------------------------------------------------------------------  
    PLATE_THICK = 25.4/8*1e-3; % 1/8 inch thick
    COIL_PLATE_SEP_DIST = 0;
    SEP_DIST = (1.0825+0.7095)*25.4e-3;
    
    % Define region
    L.layerN = 5;

    % Layer 1 (approaching from -inf)
    L.names{1} = 'Air';
    L.mu_r(1) = 1;
    L.zN(1) = 800;
    L.bnds(1) = -0.5; % boundary below
    L.sig(1) = 0;
    
    % Layer 2 (approaching from -inf)
    L.names{2} = 'Copper';
    L.mu_r(2) = 1;
    L.zN(2) = 500;
    L.sig(2) = 5.8e7;
    L.bnds(2) = COIL_PLATE_SEP_DIST-SEP_DIST-PLATE_THICK; % boundary below
    
    % Layer 3
    L.names{3} = 'Air containing coil';
    L.mu_r(3) = 1;
    L.zN(3) = 1000;
    L.bnds(3) = COIL_PLATE_SEP_DIST-SEP_DIST; % boundary below
    L.sig(3) = 0;

    % Layer 4
    L.names{4} = 'Copper';
    L.mu_r(4) = 1;
    L.zN(4) = 500;
    L.sig(4) = 5.8e7;
    L.bnds(4) = COIL_PLATE_SEP_DIST; 
    
    % Layer 5
    L.names{5} = 'Air';
    L.mu_r(5) = 1;
    L.zN(5) = 800;
    L.bnds(5) = COIL_PLATE_SEP_DIST+PLATE_THICK; % boundary below
    L.sig(5) = 0;
    
    % Upper bound
    L.bnds(6) = 0.5;

    % Excitation
    L.coil_layer = 3;
end

end


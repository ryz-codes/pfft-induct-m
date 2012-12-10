function [zl dz zu lid] = init_steps(L,w)
    %INIT_STEPS Initializes the z divisions.
    % The logic is as follows, in the stated priority chain
    % 1. Maximum size of mesh region above material should be 0.5
    % 2. Maximum number of skin depths meshed should be 10 over all layers
    % 3. Maximum element size is 0.25mm
    % 4. Pick element size to be 1/50th of skin depths.
    %
    % Outputs:
    % zl, dz, zu, 
    
end


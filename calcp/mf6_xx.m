function [ Mx ] = mf6_xx(O, L, W, H, g_xx)
    %MF6 Calculates the mutual inductance due to the x->x and y->y Green's 
    % function between two rectangular blocks using quadrature.
    %
    % Method is vectorized. This means we can pass in many many rows of
    % combinations, and the method will process them together to save CPU time.
    %
    % Inputs for p combinations of filaments:
    %    O - origin vector, p x 6 column vector
    %    L - length vector, p x 6 column vector
    %    W - width vector, p x 6 column vector
    %    H - height vector, p x 6 column vector
    %
    % The source filament occupies the first 3 columns. The field filament
    % occupies the last 3 columns. Each 3-columns are ordered [x, y, z]

    % details flag
    DETAILS = true;
    if nargin == 4
        g_xx = @green_xx;
    end
    
    %-----
    % Perform volume integration over both blocks. 
    %-----
    % Quadrature nodes and weights. Use persistent keyword to make them static
    persistent lu1
    persistent lu2
    persistent lu3
    persistent lu4
    persistent lu5
    persistent lu6
    persistent w
    if isempty(lu1) % Fetch the nodes to call once
        [lu w] = nwspgr('KPU', 6, 3);
        lu1 = lu(:,1).';
        lu2 = lu(:,2).';
        lu3 = lu(:,3).';
        lu4 = lu(:,4).';
        lu5 = lu(:,5).';
        lu6 = lu(:,6).';
        if DETAILS
        fprintf('Sparse grid quadrature KPU loaded \n\t with %g nodes per integration.\n', ...
            length(lu1));
        end
    end
    
    % timer
    if DETAILS, fprintf('===BEGIN mf6_xx===\n'); tic; end
    
    % define r vectors for G(r,rp) evaluation
    temp1 = ones(size(lu1));
    
    % source quadrature points (rp)
    rx1 = O(:,1)*temp1 + L(:,1)*lu1 + W(:,1)*lu2 + H(:,1)*lu3;
    ry1 = O(:,2)*temp1 + L(:,2)*lu1 + W(:,2)*lu2 + H(:,2)*lu3;
    rz1 = O(:,3)*temp1 + L(:,3)*lu1 + W(:,3)*lu2 + H(:,3)*lu3;
    
    % field quadrature points (r)
    rx2 = O(:,4)*temp1 + L(:,4)*lu4 + W(:,4)*lu5 + H(:,4)*lu6;
    ry2 = O(:,5)*temp1 + L(:,5)*lu4 + W(:,5)*lu5 + H(:,5)*lu6;
    rz2 = O(:,6)*temp1 + L(:,6)*lu4 + W(:,6)*lu5 + H(:,6)*lu6;

    r = sqrt((rx1-rx2).^2+(ry1-ry2).^2); % radial distances
    
    % timer
    if DETAILS, fprintf('Make columns time: %g\n', toc); tic; end
            
    % Lookup table
    integrand = g_xx(r,rz2,rz1);
    
    % timer
    if DETAILS, fprintf('Lookup table time: %g (%g pnts)\n', ...
            toc, numel(integrand)); tic; end
    
    Mx = (integrand)*w; % Complete quadrature by weighing with weights
    
    %-----
    % Adapt for the mutual inductance formula
    %    M = \iint Gxx \hat{x} \cdot dL1 \cdot dL2
    %-----
    % Get lengths
	len1 = sqrt(sum(L(:,(1:3)).^2,2));
    len2 = sqrt(sum(L(:,(4:6)).^2,2));
    
    % Obtain the x directional unit vector in fil1's frame of reference
    % for the dot product with Gxx \hat{x}
    ux = L(:,1:2);
    ux = ux./(sqrt(ux(:,1).^2+ux(:,2).^2)*[1 1]);

    % Rise angle calculations. If both are flat on the xy axis then both 
    % will be 1.
    cos_t1 = (L(:,1) .* ux(:,1) + L(:,2) .* ux(:,2))./len1; % yes
    cos_t2 = (L(:,4) .* ux(:,1) + L(:,5) .* ux(:,2))./len2; % yes
    
    % Multiply by rise angle to consolidate 
    Mx = Mx.*cos_t1.*cos_t2.*len1.*len2; 
    
    % timer
    if DETAILS, fprintf('Output stage time: %g\n', toc); end
            
end
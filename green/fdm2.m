function fdm2
% FDM uses non-uniformly spaced grid points to perform calculations.

% PARAMETERS
% Number of points per direction
N_r = 3;
N_z1 = 3;
N_z2 = 5;


r = slicespace(0,1,N_r)
z1 = linspace(0,1,N_z1)
z2 = slicespace(0,-1,N_z2)

% Assemble the R-variation blocks
% Equation implemented: d/dr x + r * d2/dr2 x = 0
% domain
indbar = (1:N_r)';
ind = [indbar-1 indbar indbar+1];
coeff = zeros(size(ind));
for ii = 2:(N_r-1) 
    coeff(ii,:) = fdcoeffF(1,r(ii),r(ind(ii,:)));
    coeff(ii,:) = coeff(ii,:) + r(ii)*fdcoeffF(2,r(ii),r(ind(ii,:)));
end

% left boundary;
ind(1,:) = ind(2,:); 
coeff(1,:) = fdcoeffF(1,r(1),r(ind(1,:))); % Neuman

% right boundary;
ind(end,:) = ind(end-1,:); 
coeff(end,:) = fdcoeffF(0,r(end),r(ind(end,:))); % Dirichlet

AR = sparse(repmat(indbar,1,3),ind,coeff);

% Assemble the Z-variation blocks
% Equation implemented: d2/dz2 = 0
indbar = (1:N_z1)';
ind = [indbar-1 indbar indbar+1];
coeff = zeros(size(ind));
for ii = 2:(N_z1-1) 
    coeff(ii,:) = fdcoeffF(2,z1(ii),z1(ind(ii,:)));
end

% top boundary;
ind(1,:) = ind(2,:);
coeff(1,1:2) = coeff(2,2:3); % Implicit Dirichlet

% bottom boundary;
ind(end,:) = ind(end-1,:); 
coeff(end,end-1:end) = coeff(end-1,end-2:end-1); % Implicit Dirichlet

AZ1 = sparse(repmat(indbar,1,3),ind,coeff);

% Assemble combined fdm system
A = kron(speye(N_z1),AR);
full(A)
doms = noends(ones(size(r')));
A = A + kron(AZ1,spdiags(doms,0,N_r,N_r));
spy(A)
full(A)

    function c = slicespace(a,b,N)
    % Forms a set of numbers from min(a,b) to max(a,b), exponentially
    % divided according to "split by half" rules, bunched together at the
    % "a" side. 
    % 
    % Examples:
    % >> slicespace(0,1,4)
    % ans =    0    0.2500    0.5000    1.0000
    % >> slicespace(1,0,4)
    % ans =    0    0.5000    0.7500    1.0000
    %       
        N = N-2;
        c = horzcat(0,2.^(0:N)/2^N);
        if a < b
            c = a + c*(b-a);
        elseif b < a
            c = 1-c(end:-1:1);
            c = b + c*(a-b);
        end
    end
    
    function a = noends(a)
        a(1) = 0;
        a(end) = 0;
    end
end
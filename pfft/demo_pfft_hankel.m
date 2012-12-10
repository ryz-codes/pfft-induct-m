clear;

%
%==========================================================================
% SIMPLE POINT-BASED N-D PFFT EXAMPLE
%==========================================================================
DIMS = 3;

NUM_ELE = 1e2; % Number of elements simulated
GRID_ELE = ones(1,DIMS) * 2^6; % Numer of grid points in each dimension
OUTER_BOUNDS = [ones(1,DIMS)*10; ones(1,DIMS)*0]; % Plus or minus bounds where elements are contained
INTERP_ORD = 2;
DIRECT_STEN = 4;
NUM_MV = 1; % Numer of times to repeat Matrix Vector product

CMP_CALCP = true;

% Toeplitz function (r-r')
green = @inv_r; % Use supplied 1/r green's function

% Hankel function (r+r')
hankel_fn = @inv_r; % Hankel part is also 1/r

% Generate centroids and charges within a box
X = (rand(NUM_ELE,DIMS))*8+1;
Q = (ones(NUM_ELE,1)-0.5);

calcp_T = @(from,to) inv_r(num2cell(X(from,:)-X(to,:),1));
calcp_H = @(from,to) inv_r(num2cell(X(from,:)+X(to,:),1));

%--------------------------------------------------------------------------
fprintf('\n\nBeginning PFFT for %g elements and %g grid points\n\n',...
    NUM_ELE,GRID_ELE(1).^DIMS);
    tot_t=tic;
p = pfft(GRID_ELE,OUTER_BOUNDS,INTERP_ORD,DIRECT_STEN);
p.init_kernel(green,hankel_fn,ones(1,DIMS));
p.init_geometry(X);
p.init_precorrect(calcp_T);
[V_pfft V_fft] = p.fastmv(Q);
tot_t = toc(tot_t);
V2 = Q;
tic
for repeater = 1:NUM_MV
    V2 = p.fastmv(V2);
end
    fprintf('-------\nMV time (%d run av): \t%g\n',NUM_MV,toc/NUM_MV); 
    fprintf('-------\nTotal time: \t\t%g\n',tot_t); 
%--------------------------------------------------------------------------

tic
if CMP_CALCP && NUM_ELE < 501 % Freeze-prevention   
    
    P_T = zeros(NUM_ELE);
    P_H = zeros(NUM_ELE);
    for ii = 1:NUM_ELE
        for ij = 1:NUM_ELE
            P_T(ii,ij) = calcp_T(ii,ij)+calcp_H(ii,ij);
            %P_H(ii,ij) = calcp_H(ii,ij);
        end
    end
    V_calcp = P_T*Q;
    fprintf('\n\nO(N^2) solve time: %g\n',toc);
    %disp([V_calcp, V_pfft])
    %disp((abs(V_calcp - V_pfft)./abs(V_calcp)));
    fprintf('\nFFT conv err:\t%g\n',mean(abs(V_calcp - V_fft)./abs(V_calcp)));
    fprintf('Pfft err:\t%g\n',mean(abs(V_calcp - V_pfft)./abs(V_calcp)));
    fprintf('Max pfft err:\t%g\n',max(abs(V_calcp - V_pfft)./abs(V_calcp)));
end
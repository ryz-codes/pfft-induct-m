clear;
load ferrite_iron.mat;
%
%==========================================================================
% SIMPLE POINT-BASED N-D PFFT EXAMPLE
%==========================================================================
DIMS = 3;

NUM_ELE = 10; % Number of elements simulated
GRID_ELE = [ones(1,2)*2^7 8]; % Numer of grid points in each dimension
OUTER_BOUNDS = [[-1 1; -1 1]*0.12; [-0.015 -0.001]]; % Plus or minus bounds where elements are contained
%OUTER_BOUNDS = [ones(1,DIMS)*10; ones(1,DIMS)*0]; % Plus or minus bounds where elements are contained
INTERP_ORD = 1;
DIRECT_STEN = 2;
NUM_MV = 1; % Numer of times to repeat Matrix Vector product

CMP_CALCP = true;

% Toeplitz function (r-r')
t_fun = @(in) tph_wrap(g,'T',in);
h_fun = @(in) tph_wrap(g,'H',in);

% Generate centroids and charges within a box
X = (rand(NUM_ELE,DIMS))*5e-3-1e-2;
Q = (ones(NUM_ELE,1)-0.5);

calcp_T = @(from,to) t_fun(num2cell(X(from,:)-X(to,:),1));
calcp_H = @(from,to) h_fun(num2cell([X(from,1:2)-X(to,1:2), X(from,3)+X(to,3)],1));

%--------------------------------------------------------------------------
fprintf('\n\nBeginning PFFT for %g elements and %g grid points\n\n',...
    NUM_ELE,GRID_ELE(1).^DIMS);
    tot_t=tic;
p = pfft(GRID_ELE,OUTER_BOUNDS,INTERP_ORD,DIRECT_STEN);
%p.init_kernel(t_fun);
%p.init_kernel(t_fun,h_fun,[0 0 1]);
p.init_kernel([],h_fun,[0 0 1]);
p.init_geometry(X);
%p.init_precorrect(calcp_T);
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
            %P_T(ii,ij) = calcp_T(ii,ij);%+calcp_H(ii,ij);
            P_H(ii,ij) = calcp_H(ii,ij);
        end
    end
    V_calcp = P_H*Q;
    fprintf('\n\nO(N^2) solve time: %g\n',toc);
    disp([V_calcp, V_pfft])
    %disp((abs(V_calcp - V_pfft)./abs(V_calcp)));
    fprintf('\nFFT conv err:\t%g\n',mean(abs(V_calcp - V_fft)./abs(V_calcp)));
    fprintf('Pfft err:\t%g\n',mean(abs(V_calcp - V_pfft)./abs(V_calcp)));
    fprintf('Max pfft err:\t%g\n',max(abs(V_calcp - V_pfft)./abs(V_calcp)));
end
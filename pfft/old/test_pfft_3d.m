clear;

%
% Known bugs: Translation of relative coordinates is buggy. Direct stencil 
% will sometimes screw up and contain the wrong coordinates, due to a poor
% translation from subscripts to indices.
%
% Translation of big grid indices to little grid indices can also be buggy. 
%
%==========================================================================
% SIMPLE POINT-BASED 3D PFFT EXAMPLE
%==========================================================================
DIMS = 3;

green = @inv_r;

NUM_ELE = 500; % Number of elements simulated
GRID_ELE = ones(1,DIMS) * 2^6;
OUTER_BOUNDS = ones(1,DIMS)*10; % Plus or minus bounds where elements are contained
INTERP_ORD = 1;
DIRECT_STEN = 2;
CMP_CALCP = true;

% Generate centroids and charges
X = (rand(NUM_ELE,DIMS)-0.5)*9;
Q = (rand(NUM_ELE,1)-0.5);

calcp = @(from,to) 1./sqrt(sum((X(from,:)-X(to,:)).^2,2)+1e-6);

%--------------------------------------------------------------------------
fprintf('\n\nBeginning PFFT for %d elements and %g grid points\n\n',...
    NUM_ELE,GRID_ELE(1).^DIMS);
    tot_t=tic;
    tic
p = pfft(GRID_ELE,OUTER_BOUNDS,INTERP_ORD,DIRECT_STEN);
    fprintf('Grid preparation: \t%g\n',toc); tic
p.init_kernel(green);
    fprintf('Kernel preparation: \t%g\n',toc); tic
p.init_geometry(X);
    fprintf('Geometry preparation: \t%g\n',toc); tic
p.init_calcp(calcp);
    fprintf('Direct evaluation: \t%g\n',toc); tic
[V_pfft V_fft] = p.fastmv(Q);
    fprintf('MV time: \t\t%g\n',toc); 
    fprintf('-------\nTotal time: \t\t%g\n',toc(tot_t)); 
%--------------------------------------------------------------------------

V_dir = p.Dmat2*Q;

tic
if CMP_CALCP && NUM_ELE < 501 % Freeze-prevention   
    P2 = zeros(NUM_ELE);
    for ii = 1:NUM_ELE
        for ij = 1:NUM_ELE
            P2(ii,ij) = calcp(ii,ij);
        end
    end
    V_calcp = P2*Q;
fprintf('\n\nO(N^2) solve time: %g\n',toc);
    % Displace
    %disp([V_calcp, V_pfft]);
    %disp([V_calcp-V_pfft]);
    fprintf('\nFFT conv err:\t%g\n',mean(abs(V_calcp - V_fft)./abs(V_calcp)));
    fprintf('Direct err:\t%g\n',mean(abs(V_calcp - V_dir)./abs(V_calcp)));
    fprintf('Pfft err:\t%g\n',mean(abs(V_calcp - V_pfft)./abs(V_calcp)));
end
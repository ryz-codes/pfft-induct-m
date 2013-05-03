function [Zt,res] = fasthenry(geom,f,xybnd,zbnd)
%FASTHENRY Free-space terminal impedance extraction akin to FastHenry. 
% The input GEOM is a struct with five fields: O, L, W and L are geometry
% constructs. GEOM.type is a string that determines the type of extraction
% performed.
%
% For a simple single terminal extraction, GEOM.type="1term". Here, O, L, 
% W and H are one-level cell structures; the lowest level are filaments, 
% bundled into segments at the upper level.
%
% For a multi-terminal, single excitation extraction, GEOM.type="xterm1".
% Here, O,L,W,H are two-level cell structures; the lowest level are
% filaments, the mid levels are segments and the upper level are the
% terminals.

% default grid bounds
if nargin < 3
    xybnd = 0.16; 
end
if nargin <4
    zbnd = [-0.015 0.005];
end

% SELECT PRECONDITIONER TYPE. 1 for diagonal L, 2 for block L.
PRECOND_ = 1;

% List of frequencies
w = 2*pi*f;

% path is either a matrix or a cell array. 
% rad is the radius of the path. 
O = geom.O; L = geom.L; W = geom.W; H = geom.H;

% Check the problem type - multiterminal or single terminal
% TODO: Add more comprehensive error checks
if (iscell(O) && iscell(L) && iscell(W) && iscell(H))
    if (iscell(O{1}) && iscell(L{1}) && iscell(W{1}) && iscell(H{1}))
        TYPE_ = 2;
    else
        TYPE_ = 1;
    end
else
    TYPE_ = 0;
end

% Check if we have a tph problem. 
TPH_ = 0; % Toeplitz plus hankel flag
if isfield(geom,'g')
    % Do error checks
    assert(iscell(geom.g),'geom.g must be a cell');
    assert(numel(geom.g)==numel(f),'geom.g must contain as many elements as there are frequencies');
    assert(isa(geom.g{1},'tph'),'Each element of geom.g must be a tph object');
    TPH_ = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expand out to branches, generate mesh matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch TYPE_ % 0 no terminal, 1 single terminal, 2 multi terminal
%--------------------------------------------------------------------------    
    case 0 % No terminal
        error('No-terminal extraction not implemented');
%--------------------------------------------------------------------------
    case 1 % 1 terminal
        % Make segment-based mesh matrix
        Nfils = size(O{1},1); % number of fils per segment
        Nsegs = numel(O);
        Nbranch = Nfils*Nsegs;
        M = meshMat(Nsegs,Nfils);
        Nmesh = size(M,1);

        % Expand out the branches
        O = vertcat(O{:}); L = vertcat(L{:}); W = vertcat(W{:}); H = vertcat(H{:});
        
        % Terminal matrix
        S = 1; %one terminal

%--------------------------------------------------------------------------
    case 2 % multi-terminal
        Nterms = numel(O);
        M = cell(Nterms,1); % Mesh matrix
        S = zeros(Nterms,1); % Terminal matrix
        S(1) = 1;
        for ij = 1:Nterms
            Nfils = size(O{ij}{1},1);
            Nsegs = numel(O{ij});
            M{ij} = meshMat(Nsegs,Nfils);
            if ij < Nterms
                % Find where the terminal is
                S(ij+1) = S(ij) + size(M{ij},1);
            end
            
            % Expand out the branches
            O{ij} = vertcat(O{ij}{:}); L{ij} = vertcat(L{ij}{:}); 
            W{ij} = vertcat(W{ij}{:}); H{ij} = vertcat(H{ij}{:});
        end
        
        % Expand out the terminals
        M = blkdiag(M{:});
        O = vertcat(O{:}); L = vertcat(L{:}); 
        W = vertcat(W{:}); H = vertcat(H{:});
        
        % Total numbers
        Nbranch = size(O,1);
        Nmesh = size(M,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store results (if they are asked for)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 2
    res.M = M; % store mesh matrix
    res.S = S; % store terminal matrix
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impedance extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display
fprintf('Beginning computation for %d elements...\n',Nbranch);

% Diagonal Resistance matrix
Rself = resist(ones(Nbranch,1),[],L,W,H);

% Set up the grid
induct(xybnd,zbnd); % Set up perturbation component
[Lself, dat2] = induct([],O,L,W,H); % Precorrect, get preconditioner
res.proj = dat2.proj;
res.prec = dat2.prec;
if PRECOND_ == 1 % Diagonal-L preconditioner
    % Isolate the main diagonal, leave as dense matrix
    Lself = spdiags(Lself,0);
elseif PRECOND_ == 2
    % Turn into a sparse matrix instantly
    Rselfsp = spdiags(Rself,0,Nbranch,Nbranch);
end

% Objective vector
Vm0 = zeros(Nmesh,1); Vm0(S) = 1;

% Solve for each frequency
Zt = zeros(size(w));
if nargout == 2 % if additional info is desired
    res.Im0 = zeros(Nmesh,numel(w)); % store mesh currents
    res.solve_time = zeros(1,numel(w));
    res.iters = zeros(1,numel(w));
    res.residual = zeros(1,numel(w));
    if TPH_ == 1
        res.kernel_time = zeros(1,numel(w));
    end
    
end
for ii = 1:length(w)
    fprintf('\nStartin f = %e\n',f(ii));
    
    % Initiate the new kernel
    if TPH_ == 1
        fprintf('Updating new kernel... '); tt = tic;
        induct(geom.g{ii});
        tt = toc(tt);
        fprintf('Done in %g\n',tt);
    end
    
    % Set up preconditioner
    % P=M(R+jwLself)MT
    if PRECOND_ == 1 % diagonal L
        % Form inverse directly
        Yself = 1./(Rself+1j*w(ii).*Lself);
        Pinv = M*spdiags(Yself,0,Nbranch,Nbranch)*M';
        
        % Solve, precorrect inside the loop
        st = tic;
        [Im0, ~,residual,it] = gmres(@MZMT,Vm0,[],1e-6,200);
        st = toc(st);
        Im0 = Pinv*Im0; % Must precorrect again.
        Zt(ii)=1/sum(Im0(S));
    elseif PRECOND_ == 2 % block L
        % Form inverse by incomplete LU factorization
        setup.droptol = 1e-1;
        P = M*(Rselfsp+1j*w(ii).*Lself)*M';
        [L, U] = ilu(P,setup);
    
        % Solve
        st = tic;
        [Im0, ~,residual,it] = gmres(@MZMT,Vm0,[],1e-6,200,L,U);
        st = toc(st);
        Zt(ii)=1/sum(Im0(S));
    end
    
    if nargout == 2
        res.Im0(:,ii) = Im0; 
        res.solve_time(ii) = st;
        res.iters(ii) = max(it);
        res.residual(ii) = residual(1);
        if TPH_ == 1
            res.kernel_time(ii) = tt;
        end
    end
    
    fprintf('\nIterations: %g, Residual: %e, Time: %d\n',max(it),residual(1),round(st));
 
end

    % Multiply by mesh impedance matrix
    function Vm = MZMT(Im)
        if PRECOND_ == 1
            Ib = M'*(Pinv*Im);
        else
            Ib = M'*Im; 
        end
        Vm = M*(1j*w(ii)*induct(Ib)+Rself.*Ib);
        fprintf('*');
    end

    % Generate the mesh matrix for impedance extraction.
    % Nsegs = number of segments. Nfils = number of filaments per segment.
    function M_ = meshMat(Nsegs_,Nfils_)
        M1 = sparse(1,(0:Nsegs_-1)*Nfils_+1,1,1,Nfils_*Nsegs_);
        M2 = [ones(Nfils_-1,1), spdiags(-ones(Nfils_-1,1),0,Nfils_-1,Nfils_-1)];
        M3 = cell(1,Nfils_);
        for i_ = 1:Nsegs_
            M3{i_} = M2;
        end
        M_ = [M1; blkdiag(M3{:})];
    end
end
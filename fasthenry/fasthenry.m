function [Zt] = fasthenry(geom,f)
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

% List of frequencies
w = 2*pi*f;

% path is either a matrix or a cell array. 
% rad is the radius of the path. 
O = geom.O; L = geom.L; W = geom.W; H = geom.H;

% Check the problem type.
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

switch TYPE_ % 0 no terminal, 1 single terminal, 2 multi terminal
%--------------------------------------------------------------------------    
    case 0 % No terminal
        error('No-terminal extraction not implemented');
%--------------------------------------------------------------------------
    case 1 % 1 terminal
        % Make segment-based mesh matrix
        Nfils = size(O{1},1); % number of fils per segment
        Nsegs = numel(O);
        M = meshMat(Nsegs,Nfils);
        Nmesh = size(M,1);

        % Expand out the branches
        O = vertcat(O{:}); L = vertcat(L{:}); W = vertcat(W{:}); H = vertcat(H{:});

        % Display
        fprintf('Beginning computation for %d elements...\n',Nfils*Nsegs);
        
        % Set up the grid
        induct(O,O+L+W+H);
        induct(ones(Nfils*Nsegs,1),O,L,W,H);
        
        % Objective vector
        Vm0 = zeros(Nmesh,1); Vm0(1) = 1;

        % Solve for each frequency
        Zt = zeros(size(w));
        for ii = 1:length(w)
            Im0 = gmres(@MZMT,Vm0,[],1e-3,50);
            Zt(ii)=1/Im0(1);
        end
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
        
        % Display
        fprintf('Beginning computation for %d elements...\n',Nbranch);
        
        % Set up the grid
        induct(O,O+L+W+H);
        induct(ones(Nbranch,1),O,L,W,H);

        % Objective vector
        Vm0 = zeros(Nmesh,1); Vm0(S) = 1;

        % Solve for each frequency
        Zt = zeros(size(w));
        for ii = 1:length(w)
            Im0 = gmres(@MZMT,Vm0,[],1e-3,50);
            Zt(ii)=1/sum(Im0(S));
        end
end

    % Multiply by mesh impedance matrix
    function Vm = MZMT(Im)
        Ib = M'*Im;
        Vm = M*(1j*w(ii)*induct(Ib)+resist(Ib,[],L,W,H));
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
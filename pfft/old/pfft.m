classdef pfft < handle
    %PFFT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Grids
        pfft_grid;
        pc_gc;
        pc_gc_sten;
        
        % Stencils
        P_sten;
        I_sten;
        pc_I_sten;
        D_sten;
        
        IP_same; % not implemented..yet
        
        % Precorrection terms
        pc_grid_conv; % Little grid for direct interactions
        pc_gid_big; pc_gid_small; % maps big grid id with little grid
        pc_self_conv; % self term convolution matrix
        pc_nei_c; 
        pc_nei_r;
        pc_nei_e;
        
        % Kernel dependent terms
        Tker; % Toeplitz kernel 
        Hker; % Hankel kernel (not implemented)
        
        % Matrices
        Pmat;
        Imat;
        Dmat;
        Dmat2;
        
        ready=0;
    end
    
methods
    function hObj = pfft(grid_ele, outerbnds, interp_ord, direct_sten)

        % Make Big grid
        d = 2*outerbnds(:)./(grid_ele(:));
        g = grid(d, grid_ele);

        % Make stencils
        P_sten = stencil(interp_ord,g); % Interp / Proj stencil
        D_sten = stencil(direct_sten,g,'round'); % Direct interaction stencil

        %%% Precorrection setup
        % Make a little grid for precorrecting the close-by interactions
        Nc = (direct_sten+interp_ord)*2+1;
        gc = grid(g.d,(Nc*ones(size(g.N)))); 
        % stencil version of our little grid. Used to do relative
        % conversions
        gc_sten = stencil(direct_sten+interp_ord,g,'square'); 

        % Make an interpolation stencil in the little grid, for indexing in
        % the new environment
        precor_sten = stencil(interp_ord,gc);

        % Save the created variables
        % Grids
        hObj.pfft_grid = g; hObj.pc_gc = gc; hObj.pc_gc_sten = gc_sten;
        % Stencils
        hObj.pc_I_sten = precor_sten; hObj.P_sten = P_sten;
        hObj.I_sten = P_sten; hObj.D_sten = D_sten;

        hObj.ready(1) = 1;
    end
        
    
    % TODO: Specify limitations on the green's function
    function init_kernel(hO, green)
        % Import variables
        g = hO.pfft_grid;
        
        % Convolution setup
        green_fn = green(g.getCoordBox); 
        
        % save the kernel in the transformed form
        hO.Tker = real(fftn(green_fn));
        
        % Make precorrection convolution matrix
        hO.pc_grid_conv = pfft.convmat(green_fn, ...
            hO.P_sten.sten, hO.pc_gc_sten.sten);
        hO.pc_self_conv = pfft.convmat(green_fn, ...
            hO.P_sten.sten, hO.I_sten.sten);
    end
    
    % Inputs:
    % Centers, quadrature points, weights, individual excitations
    function init_geometry(hO, centroid, qpnt, varargin)
        % Import variables
        g = hO.pfft_grid;
        
        % Allocate and center the points to their respective stencils
        [cen_gid cen_coeff] = g.allocate(hO.P_sten,centroid); % Og is the grid point index of each
        if nargin >2
            qpnt_centered = cell(g.dims,1);
            for ii = 1:g.dims
                qpnt_centered{ii} = bsxfun(@minus,qpnt{ii},cen_coeff(:,ii));
            end
            % Make projection matrix
            this_p = pmat(hO.P_sten,qpnt_centered,varargin{:});
        else % No quadrature points. Point-based example
            cen_centered = centroid - cen_coeff;
            % Make projection matrix
            this_p = pmat(hO.P_sten,num2cell(cen_centered,1));
        end

        
        hO.Pmat = pfft.getSparseMats(this_p, cen_gid, hO.P_sten);

        % Calculate Precorrection
        [nei_from, nei_to, nei_rel_gid, num_nei] = ... 
            pfft.findNeighbors(hO.D_sten,cen_gid); % Find neighbors

        nei_r = zeros(num_nei,1);
        nei_c = zeros(num_nei,1);
        nei_e = zeros(num_nei,1);
        m_idx = false(num_nei,1);
        pntr = 1;
        
        % Little grid terms
        gid_map = hO.pc_gc_sten.sten_gid;
        Isten_gid = hO.pc_I_sten.sten_gid;
        
        for ii = 1:numel(nei_from)  
            % Extract info
            from_id = nei_from{ii};     to_id = nei_to{ii};
            num_to = numel(to_id);      num_from = numel(from_id);
            these_gid_big = nei_rel_gid{ii}; 
            
            % Project all "from" points onto the little grid, one column
            % for each "from" point
            this_fld = hO.pc_grid_conv*this_p(:,from_id);

            % gid_map is a map for all points in the little grid
            % the little grid and their corresponding relative index
            % values in the big grid. This avoids us from converting an
            % index into subscripts and back into an index in the small
            % grid.
            [these_gid ~] = find(bsxfun(@eq,these_gid_big',gid_map));
            
            % Do each "to" individually
            theseI = this_p(:,to_id).';
            this_pc = zeros(num_to,num_from);
            for ik=1:num_to
                this_pc(ik,:) = theseI(ik,:)...
                    *this_fld(these_gid(ik)+Isten_gid,:);
            end
            
            % Store the Precorrection
            pntr_end = pntr + num_from*num_to -1;
            nei_c(pntr:pntr_end) = repmat(to_id,num_from,1); % column is to
            nei_r(pntr:pntr_end) = reshape(...
                repmat(from_id',num_to,1),[],1); % row is from
            nei_e(pntr:pntr_end) = reshape(this_pc,num_from*num_to,1);
            
            % Find the mutual terms and store index
            m_idx(pntr:pntr_end) = all(...
                bsxfun(@ne,nei_c(pntr:pntr_end),from_id'),2);
            
            pntr = pntr_end+1;
        end
        
        % Replicate mutual terms symmetrically
        hO.pc_nei_c = [nei_c; nei_r(m_idx)];
        hO.pc_nei_r = [nei_r; nei_c(m_idx)];
        hO.pc_nei_e = [nei_e; nei_e(m_idx)];

    end
        
    function init_calcp(hO, calcp)
        nei_d = calcp(hO.pc_nei_c,hO.pc_nei_r);
        hO.Dmat2 = sparse(hO.pc_nei_r,hO.pc_nei_c,nei_d); % 
        
        nei_2 = nei_d - hO.pc_nei_e; % subtract precorrection
        hO.Dmat = sparse(hO.pc_nei_r,hO.pc_nei_c,nei_2);
        hO.pc_nei_c = [];
        hO.pc_nei_r = [];
        hO.pc_nei_e = [];
    end
    
    % Fast matrix vector product
    function [out out1] = fastmv(hO,I,u)
        if nargin == 2 % Scalar style
            out1 = hO.Pmat.'*fftconv(hO.Tker,hO.Pmat*I);
        elseif nargin == 3 % Vector style
            assert(size(I,1) == size(u,1),'Must provide one unit vector for each basis function');
            assert(size(u,2) == hO.p.g.dims, 'Unit vector must have as many dimensions as the PFFT')
            
            % Dot product each dimension
            out1 = 0;
            for ii = 1:size(u,2)
                if any(u(:,ii) ~= 0) 
                    out1 = out1 + u(:,ii)...
                        .*(hO.Pmat.'*fftconv(hO.Tker,hO.Pmat*(u(:,ii).*I)));
                end
            end
        end

        % Precorrect
        out = out1+hO.Dmat*I;
    end
end
%==========================================================================
% PRIVATE STATIC METHODS
%==========================================================================
methods (Static = true, Access = private)
    
    %CONVMAT
    % Products the (b x c) convolution matrix for (b x dim) input 
    % subscripts and (c x dim) input subscripts, using the dim-dimensional
    % green's function "box".
    %
    % Inputs:
    %  - green: "box" of green's function.
    %  - input_ss: input subscripts, with one column for each dimension.
    %  - output_ss: output subscripts, with one column for each dimension.
    function mat = convmat( green, input_ss, output_ss)
        % Find dimensions and size of green's function box
        dims = size(input_ss,2);
        for ii = 1:dims
            N(1,ii) = size(green,ii);
        end

        % Use subscript based notation. Iterate through each dimension and
        % "fold" in that dimension.
        idx = cell(1,dims);
        for ii = 1:dims
            idx{ii} = 1+ bsxfun(@minus, output_ss(:,ii), input_ss(:,ii)');
            
            % Find negative indices and fold
            negidx = (idx{ii} <1);
            idx{ii}(negidx) = idx{ii}(negidx) + N(ii);
        end

        % Convert to index and output.
        if numel(idx) >1
            mat = green(sub2ind(N,idx{:}));
        else
            mat = green(idx{1});
        end
    end
    
    
    %GETSPARSEMATS 
    % Takes a column-condensed projection matrix and expands it to a 
    % sparse matrix using the provided stencil grid points.
    % 
    % Inputs:
    %  - CondP: column-condensed matrix
    %  - gid:   grid indices that each column corresponds to.
    %  - s:     contains stencil object with relative grid indices along
    %           each row.
    function [P I] = getSparseMats(condP, gid, s)
        num_nnz = nnz(condP);
        num_ele = size(condP,2);
        num_grid = prod(s.g.N);
        sten_gid = s.sten_gid;
        num_sten = numel(sten_gid);

        row = bsxfun(@plus,sten_gid(:),gid(:).');
        col = repmat((1:num_ele),num_sten,1);
        P = sparse(row,col,condP,num_grid,num_ele,num_nnz);
        I = P.';
    end
    
    %FINDNEIGHBORS
    % Find neighbors of each centroid within this stencil. Each
    % centroid has a grid point allocated already.
    %
    % Note that this list is in an upper triangular form. As the
    % neighbor relationship is reciprocal, only one of the two
    % relationships is recorded.
    %
    % Inputs:
    %  - dsten: Direct Stencil object.
    %  - gid:   Grid indices of each element.
    %
    % Outputs a list for each unique grid point:
    %  - nei_from: centroids occupying that grid point
    %  - nei_to:   centroids occupying neighbor grid poinds
    %  - nei_gid:  grid number difference between neighbor and self
    %  - num_nei:  number neighbors recorded in total
    function [nei_from nei_to nei_gid_fnd num_nei] = findNeighbors(dsten, gid)
        % Make "2D-map" of grid points to centroids
        N_cen = numel(gid);
        N_grid = prod(dsten.g.N);
        cenid = 1:N_cen;
        spmap = sparse(cenid,gid,1,N_cen,N_grid,N_cen);

        % Achieve upper triangular form by taking out the "looking back"
        % part.
        this_sten = dsten.sten_gid;
        this_sten = this_sten(this_sten >= 0); 

        % Iterate through each uniquely occupied grid point
        gid_unique = unique(gid);
        N_unique = numel(gid_unique);

        % Initialize output
        nei_from = cell(N_unique,1);
        nei_to = cell(N_unique,1);
        nei_gid_fnd = cell(N_unique,1);
        num_nei = 0;
        for ii = 1:N_unique 
            % Get this grid id
            this_gid = gid_unique(ii);

            % Where neighbors should be
            nei_gid = this_sten + this_gid;
            
            % Exclude those points outside our grid limits
            nei_gid = nei_gid(nei_gid > 0 & nei_gid <= N_grid);

            % Extract centroids from this grid point.
            this_cen = find(spmap(:,this_gid));

            % Extract centroids from neighbor grid points.
            % Also find where the neighbor gids actually are
            [this_nei, this_fnd] = find(spmap(:,nei_gid));

            % store
            nei_from{ii} = this_cen;
            nei_to{ii} = this_nei;
            nei_gid_fnd{ii} = nei_gid(this_fnd)-this_gid;
            num_nei = num_nei+numel(this_cen)*numel(this_nei);
        end
    end
end
    
end


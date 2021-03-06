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
        D_sten;
        
        IP_same; % not implemented..yet
        
        % Precorrection terms
        pc_grid_conv; % Little grid for direct interactions
        pc_orig_ind; % Origin index is skewed to reduce computations
        pc_nei_c; 
        pc_nei_r;
        pc_nei_e;
        cen_gid;
        
        % Kernel dependent terms
        Tker = 0; % Toeplitz kernel 
        Hker; % Hankel kernel 
        rotdims;
        
        % Matrices
        Pmat;
        Imat;
        PCmat;
        Dmat;
        
        Treal=1;
        ready=0;
        verbose=1; % mode 0 for clean, mode 1 for verbose. 
        debug=1; % Debug mode, stores Dmat2 and gives the non-precorrected results
    end
    
methods
    function hObj = pfft(grid_ele, outerbnds, interp_ord, direct_sten, mode_in)
        
        % Verbose mode
        if nargin == 5
            hObj.verbose = mode_in;
            hObj.debug = mode_in;
        end
        if hObj.verbose == 1
            t_verb = tic;
        end
        
        % Make Big grid
        if numel(outerbnds) == size(grid_ele,2) % zero centered
            d = 2*outerbnds(:)./(grid_ele(:));
            g = grids(d, grid_ele);
        elseif numel(outerbnds) == 2*size(grid_ele,2) % origin specified
            if size(outerbnds,1) == 2 % rotate if 2xN matrix given
                outerbnds = outerbnds';
            end
            d = (max(outerbnds,[],2)-min(outerbnds,[],2))./(grid_ele(:));
            orig = (max(outerbnds,[],2)-min(outerbnds,[],2))/2+min(outerbnds,[],2);
            g = grids(d, grid_ele,orig);
        end
            

        % Make stencils
        hObj.P_sten = stencil(interp_ord,g); % Interp / Proj stencil
        hObj.I_sten = hObj.P_sten; 
        
        %%% Precorrection setup
        if direct_sten ~= 0 
            hObj.D_sten = stencil(direct_sten,g,'square'); % Direct interaction stencil

            % Make a little grid for precorrecting the close-by interactions
            Nc = (direct_sten+interp_ord)*2+1;
            hObj.pc_gc = grids(g.d,(Nc*ones(size(g.N)))); 
            hObj.pc_gc_sten = stencil(direct_sten+interp_ord,g,'square');
        else
            hObj.D_sten = [];
        end

        % Save the grid
        hObj.pfft_grid = g;

        % Verbose mode
        if hObj.verbose == 1
            fprintf('Grid preparation: \t%g\n',toc(t_verb));
        end
    end
        
    
    % TODO: Specify limitations on the green's function
    % First time this function is called, it sets up all the Green's
    % functions
    % Second time this function is called, it stacks the new kernel on top
    % of the old one.
    % if it's called without an argument, it deletes the previously
    % initiated kernel so that they can be redefined. 
    function init_kernel(hO, t_fun, h_fun, rotdims)
        % Verbose mode
        if hO.verbose == 1
            t_verb = tic;
        end
        
        % Delete the kernels stored
        if nargin == 1
            hO.Hker = [];
            hO.Tker = 0;
            return;
        end
        
        % Import variables
        g = hO.pfft_grid;
        %gc = hO.pc_gc;
        
        % Convolution setup
        % save the kernel in the transformed form
        if nargin > 2
            if isempty(hO.Hker)
                hankel_ker = h_fun(g.getCoordBox(rotdims));
                
                % remove nans
                nan_id = isnan(hankel_ker);
                if any(nan_id(:))
                    warning('nan encountered in hankel lookup. They have been set to zero, but results may be wrong.');
                    hankel_ker(nan_id) = 0;
                end
                
                hO.Hker = fftn(hankel_ker);
                hO.rotdims = rotdims;
            else
                error('Cannot stack Hankel kernels, only Toeplitz ones.');
            end
        end
        
        if ~isempty(t_fun)
            toeplitz_ker = t_fun(g.getCoordBox); 
            
            % remove nans
            nan_id = isnan(toeplitz_ker);
            if any(nan_id(:))
                warning('nan encountered in toeplitz lookup. They have been set to zero, but results may be wrong.');
                toeplitz_ker(nan_id) = 0;
            end
            
            hO.Tker = hO.Tker + (fftn(toeplitz_ker));
            if isempty(hO.pc_grid_conv) % If precorrection hasn't been defined yet
                % Make stencil for convolution
                % "Reduced stencil" reduces evaluation of the convolution matrix
                % to only the forward half, plus the few behind needed for the
                % self-term.
                %
                % Precorrection can only be performed for Toeplitz-type functions!!
                %
                lowest_gid = min([hO.P_sten.sten_gid; hO.I_sten.sten_gid]);
                reduced_sten = hO.pc_gc_sten.sten_gid;
                reduced_sten = reduced_sten(...
                    find(reduced_sten==lowest_gid,1,'first'):end); 

                hO.pc_orig_ind = find(reduced_sten==0);
                reduced_sten = cell2mat(g.ind2ss_rel(reduced_sten));

                hO.pc_grid_conv = pfft.convmat(toeplitz_ker, ...
                     hO.P_sten.sten, reduced_sten);  
                if hO.verbose == 1
                    fprintf('Precorrection function: %s\n',char(t_fun));
                end
            else
                if hO.verbose == 1
                    fprintf('Added to Tker but not precorrected: %s \n',char(t_fun));
                end
            end
        end
        
        % DO NOT DELETE BELOW, USED AS BACKUP IN CASE THE ABOVE HAS BUG
%         gc_sten = (1:prod(gc.N))';
%         gc_sten = cell2mat(gc.coord(gc_sten));
%         gc_sten = bsxfun(@rdivide,gc_sten,gc.d);
%         
%         % Make precorrection convolution matrix
%         hO.pc_grid_conv = pfft.convmat(green_fn, ...
%             hO.P_sten.sten, gc_sten);        
%         hO.pc_orig_ind = gc.orig_ind;

        % Verbose mode
        if hO.verbose == 1
            fprintf('Kernel preparation: \t%g\n',toc(t_verb));
        end
    end
    
    % Inputs:
    % Centers, quadrature points, weights, individual excitations
    % if quadrature points and weights are given, place qpnts in cells, one
    % cell for each quadrature point.
    function init_geometry(hO, centroid, varargin)
        % Verbose mode
        if hO.verbose == 1
            t_verb = tic;
        end
        
        % Import variables
        g = hO.pfft_grid;
        
        % Allocate and center the points to their respective stencils
        [hO.cen_gid cen_coeff] = g.allocate(hO.P_sten,centroid); % Og is the grid point index of each
            
        if nargin == 2

            % Make projection matrix
            hO.Pmat = pfft.pmat(hO.P_sten,cen_coeff,centroid,[],[]);
        elseif nargin == 4 || nargin == 5
            
            % Make projection matrix
            hO.Pmat = pfft.pmat(hO.P_sten,cen_coeff,varargin{:});
        else
            error('Must init geometry either with only centroids or with centroids and quadrature points');
        end
            
        if hO.verbose == 1
            fprintf('Geometry projection: \t%g\n',toc(t_verb));
        end
    end
    
    % Parfor form of our precorrection algorithm,
    % suitable for parallel processing.
    %
    % vec_pc controls whether the direct and precorrection matrices are
    % stored together or separately.
    function init_precorrect(hO, calcp, vec_pc)
        if isempty(hO.D_sten)
            warning('Direct stencil is set to zero. Exiting.');
            return;
        end
        if isempty(hO.Tker)
            warning('No Toeplitz style kernel to precorrect. Exiting.');
            return;
        end
        this_p = hO.Pmat;
        if issparse(this_p)
            warning('Must not run init_precorrect more than once!');
            return;
        end
        if nargin < 3
            vec_pc = false;
        end
        
        % Verbose mode
        if hO.verbose == 1
            h = waitbar(0,'Finding Neighbors...');
            t_verb = tic;
        end
        
        g = hO.pfft_grid;
        gc = hO.pc_gc;
        
        % Project all points onto their own little grid.
        fldmat = hO.pc_grid_conv;

        % interpolation stencil
        Isten_gid = gc.ss2ind_rel(g.ind2ss_rel(hO.I_sten.sten_gid));
        
        % Find neighbors
        [from, to, rel_gid, guniq] = ... 
            pfft.findNeighbors(hO.D_sten,hO.cen_gid);
        
        % Convert big grid indices to little grid indices
        % neighbor indices
        rel_gid = cell2mat(rel_gid);
        rel_gid = gc.ss2ind_rel(g.ind2ss_rel(rel_gid))... 
                    + hO.pc_orig_ind;
        to_start = cumsum([1; cellfun(@numel,to)]);
        rel_gid_c = cell(numel(to_start)-1,1);
        for ii = 1:(numel(to_start)-1)
            rel_gid_c{ii} = rel_gid(to_start(ii):to_start(ii+1)-1);
        end
        clear('to_start','rel_gid');
                
        % Init
        Ntot = numel(from);
        nei_c = cell(Ntot,1);
        nei_r = cell(Ntot,1);
        nei_e = cell(Ntot,1);
        m_idx = cell(Ntot,1);
        
        % Verbose mode
        if hO.verbose == 1
            waitbar(0,h,'Precorrecting in parallel...');
        end

        for ii = 1:numel(guniq)  
            % Extract info
            from_id = from{guniq(ii)};         to_id = to{guniq(ii)};            
            num_to = numel(to_id);      num_from = numel(from_id);

            % Project all "from" points onto the little grid, one column
            % for each "from" point
            this_fld = fldmat*this_p(:,from_id);

            % Collect little grid index of the "to" points
            %these_gid = rel_gid(to_start(ii):(to_start(ii+1)-1));
            these_gid = rel_gid_c{guniq(ii)};

            % Do each "from" individually
            theseI = this_p(:,to_id); % CHANGE THIS TO this_i LATER
            to_ndx = bsxfun(@plus,these_gid',Isten_gid);
            this_pc = zeros(num_to,num_from);
            for ik=1:num_from
                my_fld = this_fld(:,ik);
                try
                my_fld = my_fld(to_ndx);
                catch
                    disp(' ');
                end
                this_pc(:,ik) = sum(theseI.*my_fld,1);
            end

            % Store the Precorrection
            nei_c{ii} = repmat(to_id,num_from,1); % column is to
            nei_r{ii} = reshape(...
                repmat(from_id',num_to,1),[],1); % row is from
            nei_e{ii} = reshape(this_pc,num_from*num_to,1);


            % Find the mutual terms and store index
            m_idx{ii} = all(...
                bsxfun(@ne,nei_c{ii},from_id'),2);
        end
        
        % Clear the memory heavy variables
        clear('from','to','rel_gid_c','fldmat')

        % remove trailing zeros
        nei_c = vertcat(nei_c{:});
        nei_r = vertcat(nei_r{:});
        nei_e = vertcat(nei_e{:});
        
        % Verbose mode
        if hO.verbose == 1
            waitbar(0,h,'Calculating direct terms...');
            fprintf('Direct terms: \t%g\n',toc(t_verb));
            t_verb = tic;
        end
        
        if vec_pc
            nei_d = nei_e;
            nei_e = zeros(size(nei_e));
        end
        
        % Calculate direct elements by blocks
        num_ent = numel(nei_e);
        if num_ent < 3e6
            nei_e = calcp(nei_c,nei_r) - nei_e;
            fprintf('Computing %d direct interactions...\n',numel(nei_e));
        else
            BLKSIZ = 3e6;
            vblks = BLKSIZ*ones(floor(num_ent/BLKSIZ),1);
            vblks = vertcat(vblks,mod(num_ent,BLKSIZ));
            
            nei_c = mat2cell(nei_c,vblks);
            nei_r = mat2cell(nei_r,vblks);
            nei_e = mat2cell(nei_e,vblks);
            fprintf('Computing direct interactions as %d x 3e6 sized blocks...\n',numel(nei_e));
            for ii = 1:numel(nei_e)
                nei_e{ii} = calcp(nei_c{ii},nei_r{ii}) - nei_e{ii};
            end

            % remove trailing zeros
            nei_c = vertcat(nei_c{:});
            nei_r = vertcat(nei_r{:});
            nei_e = vertcat(nei_e{:});
        end
        
        % Replicate mutual terms symmetrically
        m_idx = logical(vertcat(m_idx{:}));
        nei_c = [nei_c; nei_r(m_idx)];
        nei_r = [nei_r; nei_c(m_idx)];
        nei_e = [nei_e; nei_e(m_idx)]; 
        if vec_pc
            nei_d = [nei_d; nei_d(m_idx)]; 
        end
        clear('m_idx');
        
        % Verbose mode
        if hO.verbose == 1
            waitbar(0,h,'Froming matrices...');
        end
        
        % Make the sparse matrices
        hO.Dmat = sparse(nei_r,nei_c,nei_e);
        if vec_pc
            hO.PCmat = sparse(nei_r,nei_c,nei_d);
        end
        hO.Pmat = pfft.getSparseMats(this_p, hO.cen_gid, hO.P_sten);

        % Verbose mode
        if hO.verbose
            close(h);
            fprintf('Direct terms: \t%g\n',toc(t_verb));
        end
    end
    
    
    function [out, out_proj] = fastmv(hO,Q,u)
        % For whatever reason Pmat wasn't sparsified
        if ~issparse(hO.Pmat)
            hO.Pmat = pfft.getSparseMats(hO.Pmat, hO.cen_gid, hO.P_sten);
        end
        
        % FFT solution
        if nargin < 3 % Scalar convolution
            out_proj = hO.Pmat.'*conv(hO,hO.Pmat*Q);
            out = out_proj;
        else % Vector convolution
            assert(size(u,2)==hO.pfft_grid.dims);
            
            % Get rid of the components that don't have anything
            u = u(:,sum(u.^2,1) ~= 0);
            Qu = bsxfun(@times,Q,u); % Q in each component
            
            % Convolve each component
            Pmat1 = hO.Pmat;
            
            % Skips some steps if the fields are not desired.
            if nargout == 1
                out_proj = zeros(size(u));
                parfor ii = 1:size(u,2) % parfor each component
                    out_proj(:,ii) = Pmat1'*conv(hO,Pmat1*Qu(:,ii));
                end
                
                if ~isempty(hO.PCmat)
                    % Precorrect
                    out = out_proj - hO.PCmat*Qu;
                else
                    out = out_proj;
                end
            elseif nargout == 2
                out_proj = zeros(size(Pmat1,1),size(u,2));
                parfor ii = 1:size(u,2) % parfor each component
                    out_proj(:,ii) = conv(hO,Pmat1*Qu(:,ii));
                end
                
                % Interpolate onto basis functions
                if ~isempty(hO.PCmat)
                    % Precorrect
                    out = Pmat1'*out_proj - hO.PCmat*Qu;
                else
                    out = Pmat1'*out_proj;
                end
            end
            
            % Combine and output
            out = sum(bsxfun(@times,out,u),2);
        end
        
        if ~isempty(hO.Dmat) % separately precorrected
            out = hO.Dmat*Q+out;
        end
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
    %  - from: centroids occupying that grid point
    %  - to:   centroids occupying neighbor grid poinds
    %  - gid_fnd:  grid number difference between neighbor and self
    %  - num_nei:  number neighbors recorded in total
    function [from, to, gid_fnd, guniq, num_nei] = findNeighbors(dsten, gid)
        % Make "2D-map" of grid points to centroids
        gN = dsten.g.N;
        N_grid = prod(gN);
        N_cen = numel(gid);
        cenid = 1:N_cen;
        spmap = sparse(cenid,gid,1,N_cen,N_grid,N_cen);
        
        % Achieve upper triangular form by taking out the "looking back"
        % part.
        this_sten = dsten.sten;

        % Initialize output
        from = cell(N_grid,1);
        to = cell(N_grid,1);
        gid_fnd = cell(N_grid,1);
        guniq = zeros(N_grid,1);
        num_nei = 0;
        
        for gi = 1:N_grid 
            % Extract centroids from this grid point.
            [this_cen ~] = find(spmap(:,gi));
            
            if ~isempty(this_cen)
                % Fine out where this point is as a subscript
                this_sub = cell(1,numel(gN));
                [this_sub{:}] = ind2sub(gN,gi);
                this_sub = [this_sub{:}]; % compress to array
                
                % Where neighbors should be
                nei_sub = bsxfun(@plus,this_sten,this_sub);
                
                % Exclude those neighbors outside of range
                nei_ok = all(bsxfun(@lt,nei_sub,gN) & nei_sub > 1,2);
                nei_sub = num2cell(nei_sub(nei_ok,:),1);
                
                % Convert into indices
                nei_gid = sub2ind(gN,nei_sub{:});
                
                % Exclude those that look back
                nei_gid = nei_gid(nei_gid >= gi);

                % Extract centroids from neighbor grid points.
                % Also find where the neighbor gids actually are
                [this_nei, this_fnd] = find(spmap(:,nei_gid));

                % store
                from{gi} = this_cen;
                to{gi} = this_nei;
                gid_fnd{gi} = nei_gid(this_fnd)-gi;
                guniq(gi) = gi;
                num_nei = num_nei+numel(this_cen)*numel(this_nei);
            end
        end
        guniq(guniq==0) = [];
    end
    
    %PMAT Generates the interpolation / projection matrix using the POLYNOMIAL
    % method.
    % 
    % Project individual points to a stencil.
    %   P = pmat(s,qpnt)
    %
    % Inputs
    %   s         projection / interpolation stencil 
    %   qpnt      cell of projected points
    %
    % Fit basis quadrature points to stencil points
    %   P = pmat(s,qpnt,wei)  
    %   P = pmat(s,qpnt,wei,qexc) 
    %
    % Inputs for Nb basis functions, each containing 5 quadrature points
    %   qpnt      cell of basis quadrature points {(Nb x 5),(Nb x 5),(Nb x 5)}
    %   qwei      vector of quadrature weights (1 x 5)
    %   qexc      vector of basis excitation (eg length or area) (Nb x 1) 
    %
    % Outputs 
    %   P         matrix of projection weights (Nb x Ns)
    %
    % Spreads the quadrature charges onto the provided grid coordinates via
    % collocation at these grid coordinates by solving the equation
    %
    % w_g * G  = w_q * Q
    % w_g = w_q * Q * inv(G)
    %
    % G is the rect matrix of the polynomial basis functions evaluated at 
    % each grid point (rows). The contributions are weight-summed 
    % according to the row vector w_g.
    %
    % Q is the square matrix of the Green's function evaluated at each
    % collocation point (same as quad points, rows), with the source
    % charges placed at the quadrature nodes (columns). The contributions are
    % weight-summed according to the row vector w_q.
    function [P] = pmat(s,cen,qpnt,qwei,qexc)
        ord = s.L * 2; % Order is equal to the length of the stencil
        sten_coord = s.getStenCoord; % retrieve coords
        Ns = size(sten_coord,1); % number of stencil points
        dims = size(cen,2);
        Nb = size(cen,1); % number of bases 
        
        
        % TODO: add input checks

        % Get inversion of G
        invG = pinv(f_r(ord,sten_coord));

        % initialize
        P = zeros(Ns,Nb);

        % Iterate and quadrature each basis function
        if (nargin > 2) && any(qwei) % yes quadrature
            assert(iscell(qpnt),'Quadrature points must be entered as a cell')
            assert(numel(qwei) == numel(qpnt), 'Must give the same number of quadrature weights as points');
            parfor iq = 1:numel(qwei) % repeat for each quadrature point
                P = P + qwei(iq)*(f_r(ord,qpnt{iq}-cen)*invG).';
            end
        else % no quadrature
            P = (f_r(ord,qpnt-cen)*invG).';
        end
        
        if (nargin > 3) && any(qexc)
            assert(numel(qexc)==Nb, 'Must excite with a vector of the same length as qpnt');
            P = bsxfun(@times,P,qexc(:).');
        end

    end
end
    
end
%F_R Polynomial interpolation function
%   Gives a row vector of the interpolation basis functions
%   pnts must be column vectors, one column for each dimension.
function [out] = f_r(ord, pnts)
    persistent exponent
    persistent prev_ord

    dims = size(pnts,2);
    num_ele = size(pnts,1); 
    if ~any(any(exponent)) || (ord ~= prev_ord) || (size(exponent,1)~=dims) % we don't have exponent
        % make exponent
        % Use a base conversion algorithm to work out the exponents of the
        % polynomial scheme for the interpolation function.
        num_exp_all = (ord+1)^dims; % Number of exponents in total
        b = ord+1; % Conversion base, number of exponents per dimension
        d = 0:(num_exp_all-1); % converted numbers
        s_all = zeros(dims, num_exp_all);
        for n = dims:-1:1
            s_all(n,:) = rem(d,b);
            d = floor(d/b);
        end

        % Figure out the sum and select only those with the same order or lower
        % than our desired order
        s_sel = sum(s_all,1)<=ord;
        exponent = [];
        for ii = 1:dims
            temp = s_all(ii,:);
            exponent(ii,:) = temp(s_sel);
        end

        % Store
        prev_ord = ord;
    end

    % Multiply and generate f_r
    num_exp = size(exponent,2);
    out = 1;
    for ii = 1:dims
        out = out .* (repmat(pnts(:,ii),1,num_exp))...
            .^(repmat(exponent(ii,:),num_ele,1));
    end
end
function field = conv(hO,Q)
    persistent idx;        
    
    if any(hO.rotdims)
        Q_sav = Q;
    end
    
    % Size of the box / cube considered
    siz = hO.pfft_grid.N;
    dims = hO.pfft_grid.dims;

    % Internal function. Charge is always a vector.
    % Reshape and transform.
    if dims == 1
        Q = fft(Q(:),2*numel(Q));
    else
        Q = fftn(reshape(Q,siz),siz*2);
    end

    % Toeplitz part
    if ~isempty(hO.Tker)
        field = Q.*hO.Tker;
    else
        field = 0;
    end

    % Convolve Hankel part
    if any(hO.rotdims)
        %Q = conj(Q); % rotate everything into the opposite quadrant

        % Rotate the Q
        idx = repmat({':'}, 1, dims);
        for this_dim = 1:dims
            if hO.rotdims(this_dim)
                idx{this_dim} = siz(this_dim):-1:1;
            end
        end
        if dims == 1
            Q = fft(Q_sav(idx{:}));
        else
            Q = reshape(Q_sav,siz);
            Q = fftn(Q(idx{:}),siz*2);
        end
        % Convolve, fft is not symmetric.
        field = ifftn(field+Q.*hO.Hker,'nonsymmetric');

    elseif ~isempty(hO.Tker) % no Hankel part, fft is symmetric
        % Convolve Toeplitz part
        field = ifftn(field,'nonsymmetric');
    else
        error('Setup incomplete. No convolution performed.');
    end

    % Remove the extra rows
    if ~isempty(idx) || numel(idx) ~= dims
        idx = cell(1,dims);
        for k = 1:dims
            idx{k}   = 1:siz(k);
        end
    end
    field = (field(idx{:}));
    field = field(:);
end
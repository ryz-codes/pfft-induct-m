classdef grids < handle
    %GRID grid(d,N,O)Summary of this class goes here
    % To do: Add documentation at the top
    %   ALL INDEXING START AT 1. [1, 1, 1] is the origin!!!
    
    properties(GetAccess = 'public', SetAccess = 'private')
        d = [10 10 10] * 1e-3; % meters
        N = 2.^[3 3 3]; 
        dims = 3;
        orig = [0 0 0];
        orig_ss = {5 5 5}; % for relative conversion
        orig_ind = 293;
    end
    
    properties(GetAccess = 'private', SetAccess = 'private')
        cnr = -3.5*10*1e-3*[1 1 1]; % bottom corner where index is (0,0,0)
        
    end
    
    %PUBLIC METHODS
    methods 
%==========================================================================
% CONSTRUCTOR
%==========================================================================
    function hObj = grids(d,N,O)
        if nargin == 0, return; end
        if nargin >3
            error('grid: Only 3 parameters in constructor');
        end

        % Find the number of dimensions
        hObj.dims = numel(d);
        hObj.d = d(:).'; 

        % Other defaults
        if nargin < 3
            O = zeros(1,hObj.dims);
        end
        if nargin < 2
            N = 8 * ones(1,hObj.dims);
        end

        
        % Assert that the dimensions all match
        assert(numel(N) == hObj.dims, ...
            sprintf('Must use same dimensions\n\t numel(N)=%d, dims=%d',numel(N),hObj.dims));
        assert(numel(O) == hObj.dims, ...
            sprintf('Must use same dimensions\n\t numel(O)=%d, dims=%d',numel(N),hObj.dims));

        hObj.N = N(:).'; 
        
        hObj.cnr = O(:).' - (hObj.d .* (hObj.N-1))*0.5;
        hObj.orig = O(:).';
        
        hObj.orig_ss = hObj.coord2ss(num2cell(hObj.orig));
        hObj.orig_ind = hObj.coord2ind(num2cell(hObj.orig));
    end
        
%==========================================================================
% CONVERSION METHODS (All relative methods stop working >1 D). Must fix!
%==========================================================================
% SUBSCRIPTS <-> INDEX
    % Converts absolute subscripts to indices within this grid
    function ind = ss2ind(hObj,cin)

        siz = hObj.N;

        % Stolen from sub2ind with protection removed.
        %Compute linear indices
        k = [1 cumprod(siz(1:end-1))];
        ind = 1;
        s = size(cin{1}); %For size comparison
        for i = 1:hObj.dims,
            v = cin{i};
            %%Input checking
            if ~isequal(s,size(v))
                %Verify sizes of subscripts
                error(message('MATLAB:grid.ss2ind:SubscriptVectorSize'));
            end
            ind = ind + (v-1)*k(i);
        end
    end

    % convert absolute index into subscripts within this grid
    function ss = ind2ss(hObj,ndx)

        % Stolen from ind2sub with protection removed.
        siz = hObj.N;
        n = hObj.dims;
        k = [1 cumprod(siz(1:end-1))];
        for i = n:-1:1,
          vi = rem(ndx-1, k(i)) + 1;         
          vj = (ndx - vi)/k(i) + 1; 
          ss{i} = vj; 
          ndx = vi;     
        end
    end
    
    % converts relative absolute index into relative subscripts within this grid
    function ss = ind2ss_rel(hObj,ndx)
        % Shift everything by a quarter of a quadrant, convert like normal.
            ss = ind2ss(hObj,ndx+hObj.orig_ind);
            for ii = 1:numel(ss)
                ss{ii} = ss{ii} - hObj.orig_ss{ii};
            end
    end

    % Converts relative subscripts to relative indices within this grid!
    function ind = ss2ind_rel(hObj,cin)
        siz = hObj.N;

        % Stolen from sub2ind with protection removed.
        %Compute linear indices
        k = [1 cumprod(siz(1:end-1))];
        ind = 0;
        s = size(cin{1}); %For size comparison
        for i = 1:hObj.dims,
            v = cin{i};
            %%Input checking
            if ~isequal(s,size(v))
                %Verify sizes of subscripts
                error(message('MATLAB:grid.ss2ind:SubscriptVectorSize'));
            end
            ind = ind + (v)*k(i);
        end
%         for ii = 1:numel(cin)
%             cin{ii} = cin{ii} + hObj.orig_ss{ii};
%         end
%         ind = ss2ind(hObj,cin)-hObj.orig_ind;
    end
    
    % Converts relative subscripts to indices within this grid
    % Input is shaped as column vectors, with one column for each dimension
    %
    % Used extensively by stencil functions, who store stencil relative
    % subscripts as column lists
    function ind = ss2ind_rel_col(hObj,vin)
        siz = hObj.N;
        k = [1 cumprod(siz(1:end-1))];
        ind = vin *k.';
    end
%--------------------------------------------------------------------------
% COORDS -> INDEX / SUBSCRIPTS
    % convert coordinates to nearest subscript
    function cout = coord2ss(hObj,cin)
        % Assign and round
        cout = cell(1,hObj.dims);
        for ii=1:hObj.dims
            temp = (cin{ii}-hObj.cnr(ii)) / hObj.d(ii);
            cout{ii} = round(temp)+1;
        end
    end

    % convert coordinates to nearest index
    function ind = coord2ind(hObj,cin)
        ss = hObj.coord2ss(cin);
        ind = hObj.ss2ind(ss);
    end
    
%--------------------------------------------------------------------------
% INDEX / SUBSCRIPTS -> COORDS
    % convert index or subscript into absolute coordinates
    function cout = coord(hObj,cin)
        if iscell(cin) % subscript based
            hObj.assert_ss(cin); % Error test the input
            cout = hObj.ss2coord(false, cin);
        else % index based
            hObj.assert_ind(cin); % Error test the input
            cout = hObj.ind2coord(false, cin);
        end
    end

    % convert index or subscript into relative coordinates
    function cout = coord_rel(hObj,cin)
        if iscell(cin) % subscript based

            cout = hObj.ss2coord(true, cin);
        else % index based
            cout = hObj.ind2coord(true, cin);
        end
    end

%==========================================================================
% PFFT-SPECIFIC METHODS
%==========================================================================
    % Gives the coordinate box that can be used to generate the Green's
    % function convolution matrix
    function cout = getCoordBox(hObj,rotdims)
        this_N = hObj.N; % Get bounds
        this_d = hObj.d;
        this_dims = hObj.dims;
        if nargin == 1
            rotdims = zeros(1,this_dims); % default to not rotating anything.
        end
        assert(numel(rotdims) == this_dims,'rotdims must have the same number of dimensions as the box');
        
        % If we rotate anything then we will both with finding the corners
        if any(rotdims)
            gtop = cell2mat(hObj.coord(prod(this_N)))./this_d;
            gbot = cell2mat(hObj.coord(1))./this_d;
        end
        
        % guide is the box for this dimension
        guide = cell(this_dims,1);
        for ii = 1:this_dims
            if rotdims(ii)
                % Hankel style
                %temp = (2*gbot(ii)):1:(2*gtop(ii)+1);
                temp = (2*gbot(ii)) + (0:(2*this_N(ii)-1));
                guide{ii} = this_d(ii)*temp([this_N(ii):end, 1:(this_N(ii)-1)]).';
            else
                % Toeplitz style
                guide{ii} = this_d(ii)*[0:(this_N(ii)), (this_N(ii)-1):-1:1].';
            end
        end

        % Repeat and expand each box guide to a total box
        if this_dims>1
            M = cell(1,this_dims);
%             exec_str1 = '';
%             exec_str2 = '';
%             for ii = 1:this_dims
%                 exec_str1 = [exec_str1 sprintf('M{%d},',ii)];
%                 exec_str2 = [exec_str2 sprintf('guide{%d},',ii)];
%             end
%             exec_str = ['[' exec_str1(1:(end-1)) ... end-1 to get rid of comma
%                 '] = ndgrid(' exec_str2(1:(end-1)) ');'];
%             eval(exec_str);
            [M{:}] = ndgrid(guide{:});

            cout = M;
        else
            cout = guide;
        end
    end
    
    % Allocate a list of basis function centroids to grid indices
    function [g_ind g_coord] = allocate(hObj, s, centroid)
        % Check dimensions
        [N dims] = size(centroid);
        assert(dims == hObj.dims,'Num of cols in centroid list must be the same as the dim of the grid');

        % Check to see that all centroids are within bounds
        for ii = 1:dims
            temp1 = (centroid(:,ii) <= s.lower(ii));
            temp2 = (centroid(:,ii) >= s.upper(ii));
            if (any(temp1))                
                fprintf('Faulty centroid: %g\n',centroid(temp1,ii));
                fprintf('Lower bound: %g\n',s.lower(ii));
                error('stencil.allocate:Grid not big enough to contain all of the centroids provided');
            end
            if (any(temp2))                
                fprintf('Faulty centroid: %g\n',centroid(temp2,ii));
                fprintf('Upper bound: %g\n',s.upper(ii));
                error('stencil.allocate:Grid not big enough to contain all of the centroids provided');
            end
        end

        % convert centroid matrix to cell
        centroid_cell = mat2cell(centroid,N,ones(1,dims));

        % Allocate a grid point to each centroid
        g_ind = hObj.coord2ind(centroid_cell);

        if nargout == 2
            % Give a list of stencil coordinates
            g_coord = hObj.coord(g_ind);
            g_coord = cell2mat(g_coord);
        end
    end
    
    % Dumps all subscripts in this grid into a list
%     function mat = ss_dump(hObj)
%         this_N = hObj.N; % Get bounds
%         indx = 1:prod(this_N); % make list of the indices
%         mat = hObj.ss2ind_rel_col
%     end
    end
    
%==========================================================================
% PRIVATE METHODS
%==========================================================================
    methods (Access = private)
        
        % Convert index into coordinates
        function cout = ind2coord(hObj,rel,ind)            
            cout = hObj.ss2coord(rel,hObj.ind2ss(ind));
        end
        
        % Convert subscripts into coordinates
        function cout = ss2coord(hObj,rel,ss)
            % make cell
            cout = cell(1,hObj.dims);
            
            % for each dimension
            s = size(ss{1}); % for comparison 
            for dim = 1:hObj.dims
                v = ss{dim};
                %%Input checking
                if ~isequal(s,size(v))
                    %Verify sizes of subscripts
                    error(message('MATLAB:grid.ss2coord:SubscriptVectorSize'));
                end
                
                this_d = hObj.d(dim);
                this_cnr = 0;
                if rel == false;
                    this_cnr = hObj.cnr(dim);
                    v = v-1;
                end
                
                cout{dim} = v*this_d + this_cnr;
            end
        end
        
        % Error test a subscript
        function assert_ss(hObj,c)
            assert(numel(c) == hObj.dims, 'Number of cells must always match number of dimensions');
            siz = size(c{1});
            for ii = 1:hObj.dims
                assert(isequal(siz,size(c{ii})),'Dimensions of each cell must be the same');
                assert(min(min(c{ii}))>=1, 'Subscript out of bounds, lowest ss should be >= 1');
                assert(max(max(c{ii}))<=hObj.N(ii), 'Subscript out of bounds, largest ss should be <= N');
            end
        end
        
        % Error test an index
        function assert_ind(hObj,ind)
            assert(min(min(ind))>=1, 'Subscript out of bounds, lowest ind should be >= 1');
            assert(max(max(ind))<=prod(hObj.N), 'Subscript out of bounds, largest ind should be <= prod(N)');
        end
    end
    
%--------------------------------------------------------------------------
% STATIC TEST FUNCTIONS
%--------------------------------------------------------------------------
    methods (Static = true)
        function test
            
            fprintf('\n\nBegin grid.test():\n');
            
            % Test 4 dimensions in total
            DIMS = 4;
            D = 0.1; % space out by 0.1
            N = 8; % 8 cells in each dimension
            
            for dim = 1:DIMS
                fprintf('Testing dimension %d\n',dim);
                
                % Constructor 
                fprintf('\tConstructor...');
                g = grid(D*ones(1,dim),N*ones(1,dim));
                fprintf('OK\n');
                
                fprintf('\tTest origin is centered correctly...');
                gtop = cell2mat(g.coord(prod(g.N)));
                gbot = cell2mat(g.coord(1));
                if (isequal(gtop,-gbot) && sum(g.orig) == 0)
                    fprintf('OK\n');
                else
                    fprintf('ERROR\n');
                end
                
                
                % ind2coord
                fprintf('\tind2coord <-> coord2ind...');
                ind = 1:(N^dim);
                co = g.coord(ind);
                ind2 = g.coord2ind(co);
                if isequal(ind,ind2)
                    fprintf('OK\n');
                else
                    fprintf('ERROR\n');
                end
                
                % ind2ss
                fprintf('\tind2ss <-> ind2coord + coord2ss...');
                ind = 1:(N^dim);
                ss = g.ind2ss(ind);
                ss2 = g.coord2ss(co);
                if isequal(ss,ss2)
                    fprintf('OK\n');
                else
                    fprintf('ERROR\n');
                end
                
            end 
        end
        function test_allocate
            fprintf('\n\nBegin grid.test_allocate():\n');
            
            % Tests for stencil.allocate()
            DIMS = 4;
            SIZE = 16;
            
            fprintf('\n\nTesting grid.allocate() error handling:\n\n');
            for dim = 1:DIMS
                fprintf('Attemping dim=%d\n',dim);
                
                g = grid(ones(1,dim),SIZE*ones(1,dim));
                s = stencil(1,g);
                
                withinlims = (-(SIZE/2-2):(SIZE/2-2)).'*ones(1,dim);
                abovelims = ((SIZE/2+1):SIZE/2+4).'*ones(1,dim);
                belowlims = -abovelims;
                onupper = s.upper;
                onlower = s.lower;
            
                % This one should succeed
                fprintf('\tTrying to allocate centroid within limits. Should pass...');
                try
                    gid = g.allocate(s,withinlims);
                    fprintf('OK\n');
                catch
                    fprintf('ERROR\n');
                end

                % This one should not
                fprintf('\tTrying to allocate centroid on the limit, \n\tmust assert error...\n');
                try
                    g.allocate(s,onlower);
                    fprintf('\tERROR - NOT CAUGHT')
                catch
                    fprintf('\tOK - upper caught')
                end
                try
                    g.allocate(s,onupper);
                    fprintf('\tERROR - NOT CAUGHT\n')
                catch
                    fprintf('\tOK - lower caught\n')
                end

                % This one should not
                fprintf('\tTrying to allocate centroid with some outside the limit, \n\tmust assert error...\n');
                try
                    g.allocate(s,[belowlims;withinlims]);
                    fprintf('\tERROR - NOT CAUGHT')
                catch
                    fprintf('\tOK - upper caught')
                end
                try
                    g.allocate(s,[abovelims;withinlims]);
                    fprintf('\tERROR - NOT CAUGHT\n')
                catch
                    fprintf('\tOK - lower caught\n')
                end
            end
            
            fprintf('\n\nTesting grid.allocate() functionality:\n\n');
            for dim = 1:DIMS
                fprintf('Attemping dim=%d\n',dim);
                
                g = grid(ones(1,dim),SIZE*ones(1,dim));
                s = stencil(1,g);
                
                % pick some coords
                test_gid = s.g.coord2ind({0,0,0,0});
                test_coord = cell2mat(s.g.coord(test_gid));
                d = s.g.d; % change per cell
                
                % Test that the grid point is allocated back unto itself
                fprintf('\tTesting to see if the same point is allocated back to itself...');
                if test_gid == g.allocate(s,test_coord)
                    fprintf('OK\n');
                else
                    fprintf('ERROR\n');
                end
                
                % Test the within-limit allocation goes back to its own
                % point
                fprintf('\tTesting to see if pnts within the square is all allocated here...\n');
                for ii = 1:dim
                    moved_coord = zeros(1,dim);
                    moved_coord(ii) = 0.499*d(ii);
                    moved_coord = test_coord + moved_coord;
                    if test_gid == g.allocate(s,moved_coord)
                        fprintf('\tOK - Forward');
                    else
                        fprintf('\tERROR - Forward');
                    end
                    
                    moved_coord = zeros(1,dim);
                    moved_coord(ii) = 0.499*d(ii);
                    moved_coord = test_coord - moved_coord;
                    if test_gid == g.allocate(s,moved_coord)
                        fprintf('\tOK - Reverse dim %d\n',ii);
                    else
                        fprintf('\tERROR - Reverse dim %d\n',ii);
                    end 
                end
                
                % Test the within-limit allocation goes back to its own
                % point
                fprintf('\tTesting to see if pnts outside the square are allocated away...\n');
                for ii = 1:dim
                    moved_coord = zeros(1,dim);
                    moved_coord(ii) = 0.501*d(ii);
                    moved_coord = test_coord + moved_coord;
                    if test_gid ~= g.allocate(s,moved_coord)
                        fprintf('\tOK - Forward');
                    else
                        fprintf('\tERROR - Forward');
                    end
                    
                    moved_coord = zeros(1,dim);
                    moved_coord(ii) = 0.501*d(ii);
                    moved_coord = test_coord - moved_coord;
                    if test_gid ~= g.allocate(s,moved_coord)
                        fprintf('\tOK - Reverse dim %d\n',ii);
                    else
                        fprintf('\tERROR - Reverse dim %d\n',ii);
                    end 
                end
            end
        end
    end
end


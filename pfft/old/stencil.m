classdef stencil
    %STENCIL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        g; % Overall grid
        L; % Stencil length
        type; % Type
        sten; % Stencil in relative coordinates
        sten_gid; % grid index of sten
        
        upper; % upper limit on stencil coordinates
        lower; % lower limit on stencil coordinates
    end
    
    methods
%==========================================================================
% CONSTRUCTOR
%==========================================================================
    function hObj = stencil(L,g,type)
        if nargin == 0
            L = 1; end
        if nargin <= 1
            g = grid; end
        if nargin <= 2
            type = 'c'; end % Default to compact stencil

        % Assertions
        assert(numel(L) == 1, 'length of stencil should be scalar')
        assert(isa(g,'grid'), 'Second argument must be a grid object');

        % Figure out the type of stencil
        if strncmp(type,'s',1) || strncmp(type,'S',1) 
            hObj.sten = stencil.square(L,g.dims);
            hObj.type = 'Square';
        elseif strncmp(type,'r',1) || strncmp(type,'R',1)
            hObj.sten = stencil.round(L,g.dims);
            hObj.type = 'Round';
        elseif strncmp(type,'c',1) || strncmp(type,'C',1)
            hObj.sten = stencil.comp(L,g.dims);
            hObj.type = 'Compact';
        end

        % Convert to grid id form
        hObj.sten_gid = g.ss2ind_rel(num2cell(hObj.sten,1));

        % Find upper stencil corner
        hObj.upper = cell2mat(g.coord(prod(g.N))) - g.d*(L-0.5);

        % Find lower stencil corner
        hObj.lower = cell2mat(g.coord(1)) + g.d*(L-0.5);

        % Save stencil order
        hObj.L = L;

        % Save grid handle
        hObj.g = g;
    end

        
%==========================================================================
% PFFT-SPECIFIC METHODS
%==========================================================================
    %GETSTENCOORD
    % Gives the physical coordinates of the stencil, which can then be used
    % in projection / interpolation
    function [cout] = getStenCoord(hObj)
            
        dims = hObj.g.dims;
        temp1 = hObj.sten; % fix the silly 1 indexing issue
        Ns = size(temp1,1);

        temp1 = mat2cell(temp1,Ns,ones(1,dims));
        cout = cell2mat(hObj.g.coord_rel(temp1));
    end
        

    end
    
    
%==========================================================================
% PRIVATE METHODS
%==========================================================================
    % Static methods used within our constructor
    methods (Static = true, Access = private)
        %COMP_STEN 
        % Generates a compact stencil of a given size
        function [ cout ] = comp(L,dims)
            % First dimension
            M1 = -L:1:L;
            if dims == 1
                cout = M1.';
                return;
            end

            % Make nd grids
            M = cell(1,dims);
            exec_str = '[';
            for ii = 1:dims
                exec_str = [exec_str sprintf('M{%d} ',ii)];
            end
            exec_str = [exec_str '] = ndgrid(M1);'];
            eval(exec_str);

            % Measure the displacement distance
            d = 0;
            for ii = 1:dims
                d = d + abs(M{ii});
            end

            % Ones to keep
            keep = d<=L;

            % Output their coords
            cout = zeros(nnz(keep),dims);
            for ii = 1:dims
                temp = M{ii}(keep);
                cout(:,ii) = temp(:);
            end
        end
        
        %ROUND 
        % Generates a round stencil of a given size
        function [ cout ] = round(L,dims)
            % First dimension
            M1 = -L:1:L;
            if dims == 1
                cout = M1.';
                return;
            end

            % Make nd grids
            M = cell(1,dims);
            exec_str = '[';
            for ii = 1:dims
                exec_str = [exec_str sprintf('M{%d} ',ii)];
            end
            exec_str = [exec_str '] = ndgrid(M1);'];
            eval(exec_str);

            % Measure the radial distance
            r = 0;
            for ii = 1:dims
                r = r + M{ii}.^2;
            end
            r = sqrt(r);

            % Ones to keep
            keep = r<=L;

            % Output their coords
            cout = zeros(nnz(keep),dims);
            for ii = 1:dims
                temp = M{ii}(keep);
                cout(:,ii) = temp(:);
            end
        end
        
        %SQUARE_STEN 
        % Generates a square stencil of a given size
        function [ cout ] = square(L,dims)
            % First dimension
            M1 = -L:1:L;
            if dims == 1
                cout = M1.';
                return;
            end

            % Make nd grids
            M = cell(1,dims);
            exec_str = '[';
            for ii = 1:dims
                exec_str = [exec_str sprintf('M{%d} ',ii)];
            end
            exec_str = [exec_str '] = ndgrid(M1);'];
            eval(exec_str);

            % Output their coords
            cout = zeros(numel(M{1}),dims);
            for ii = 1:dims
                temp = M{ii};
                cout(:,ii) = temp(:);
            end
        end
    end
    
%--------------------------------------------------------------------------
% STATIC TEST FUNCTIONS
%--------------------------------------------------------------------------
    methods (Static = true)
        function test
            fprintf('\n\nBegin stencil.test():\n');
            stencil.test_allocate
            stencil.test_neighbor
        end
        
        %function test_sphere
            
        %end
        
        function test_neighbor
            fprintf('\n\nBegin stencil.test_neighbor():\n');
            
            DIMS = 4;
            for dim = 1:DIMS
            
                % Performance test
                fprintf('\n\nTesting performance:\n');

                NUM_FIL = 10000; % CHANGE THIS if needed
                NUM_GRID = 2^6; % CHANGE THIS if needed
                NUM_STEN = 3; % CHANGE THIS if needed

                fprintf('\tInitializing %d centroids, %g grid points, stencil size %d\n',...
                    NUM_FIL,NUM_GRID^3,NUM_STEN);

                g = grid(ones(1,dim),NUM_GRID*ones(1,dim));
                s = stencil(NUM_STEN,g);
                cen = (rand(NUM_FIL,dim)-0.5)*NUM_GRID*0.5;
                gid = s.allocate(cen);

                fprintf('\tRunning stencil.findNeighbors...\n');
                tic
                [nei_from nei_to] = s.findNeighbors(gid);
                fprintf('\tDONE in %g seconds\n',toc);

                % Accuracy test (to do)
                fprintf('Testing Accuracy:\n');
                NUM_FIL = 400;
                cen = (rand(NUM_FIL,dim)-0.5)*NUM_GRID*0.5;
                [gid gcoo] = s.allocate(cen);
                [nei_from nei_to] = s.findNeighbors(gid);

                % Remake neighbor matrix from nei_from and nei_to
                neigh2 = zeros(NUM_FIL);
                for ii = 1:length(nei_from)
                    from_id = nei_from{ii};
                    to_id = nei_to{ii};
                    for ij = 1:length(from_id)
                        neigh2(from_id(ij),to_id) = 1;
                        neigh2(to_id,from_id(ij)) = 1;
                    end
                end

                % Brutually find out the distance between the centroids
                DIST = zeros(NUM_FIL);
                for ii = 1:NUM_FIL
                    for ij = 1:NUM_FIL
                        DIST(ii,ij) = norm(gcoo(ii,:) - gcoo(ij,:));
                    end
                end
                % Brutually make neighbor matrix
                neigh1 = (DIST <= NUM_STEN);

                if(0 == sum(sum(neigh1-neigh2)))
                    fprintf('\tSuccessfully tested neighbor location for %d centroids\n',...
                    NUM_FIL);
                else
                    fprintf('\tFAILED testing neighbor location for %d centroids\n',...
                        NUM_FIL);
                end
            end
        end
        
        
    end
end


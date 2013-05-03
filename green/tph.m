classdef tph
%TPH - Generates a Green's function look-up table
% Generates the three-dimensional multilayered Green's function as a 
% function of r,z,zp and splits it into two dimensional Toeplitz and 
% Hankel components. This is both performed to save space and to allow fft
% techniques to be used.
%
% Syntax:  g = tph(L,lbnd,ubnd)
%
% Inputs:
%    L    - the definitions of the multilayered structure (see defaultL for 
%           descriptions of fields)
%
% Outputs:
%    g    - A green's function table object. The table can be accessed by
%           typing the following command:
%                   G = g.lookup(r,z,zp);
%
% Other m-files required: fdm_run (the fdm algorithm that generates the table)
% Subfunctions: none
% MAT-files required: none
%
% See also: defaultL.m
%
% Author: Richard Y Zhang
% Massachusetts Institute of Technology
% Email: ryz@mit.edu  
% Oct 2012; Last revision: 8th Oct 2012
    
    properties
        Hobj; Tobj; L;
        
        z;
        zp;
        r;
    end
    
    methods
        function hO = tph(L,zp,z)
            if nargin == 2
                z = zp; % centered box
            end
            % Assert that both zbnd and zpbnd have the same, even number of points.
            assert(numel(z) == numel(zp),'number of points in z must be same as zp');
            assert(mod(numel(z),2)==0,'number of points must be even');
            [H, T, r, zh, zt] = tableTPH(L,zp,z);
             
             % Make interpolant objects
             [rg1,zg1] = ndgrid(r,zh);
             [rg2,zg2] = ndgrid(r,zt);
             hO.Hobj = griddedInterpolant(rg1,zg1,H.','linear');
             hO.Tobj = griddedInterpolant(rg2,zg2,T.','linear');
             
             % Save Guides
             hO.z = z;
             hO.zp = zp;
             hO.r = r;
             
             % Layered structure
             hO.L = L;
        end
        
        function [G] = lookup(hO,r,z,zp)
            %LOOKUPTPH Summary of this function goes here
            %   Table details
            %   r,z,zp - matrix of r, z, zp of the points to be examined.
            Gh = hO.lookupH(r,z,zp);
            Gt = hO.lookupT(r,z,zp);

            G = Gh+Gt;
        end
        
        function [Gt] = lookupT(hO,r,z,zp)
            %LOOKUPTPH Summary of this function goes here
            %   Table details
            %   r,z,zp - matrix of r, z, zp of the points to be examined.
            xt = abs(z-zp);

            if numel(xt) == 1
                xt = repmat(xt,size(r));
            end
            if numel(r) == 1
                r = repmat(r,size(xt));
            end

            Gt = hO.Tobj(r,xt);
        end
        
        function [Gh] = lookupH(hO,r,z,zp)
            %LOOKUPTPH Summary of this function goes here
            %   Table details
            %   r,z,zp - matrix of r, z, zp of the points to be examined.
            xh = z+zp;

            if numel(xh) == 1
                xh = repmat(xh,size(r));
            end
            if numel(r) == 1
                r = repmat(r,size(xh));
            end
            
            Gh = hO.Hobj(r,xh);
        end
    end
    
end

function [Hr Tr r zh zt] = tableTPH(L,zp,z)
%TABLETPH Given the layer structure L, makes a T+H lookup table for the
%layer selected by the flag coil_layer.
% April 25th 2013, Added modifications based on proofs in paper. Pick
% columns 1, N, round(1.6*log(N)), c+N/2
INTERP_COLS = 4;

% Identify coil layer
coil_layer = L.coil_layer;     
zN = L.zN(coil_layer);

% Assert that this layer can contain the coil
assert(L.sig(coil_layer)==0 &&...
    L.mu_r(coil_layer)==1, 'Coil layer must be free-space!');

% Identify upper and lower boundaries
bnd_a = L.bnds(coil_layer+1);
bnd_b = L.bnds(coil_layer);

% Distribute the sources according to our proof
N = length(zp);
ip = [1,N,round(1.6*log(N))]; ip(4) = ip(3) + N/2;

tstart = tic;
A = cell(4,1);
for ii = 1:4
    fprintf('Run %u/%u: zp=%g... ',ii,4,zp(ip(ii))); tic;
    [A{ii},r] = fdm_run(z,zp(ip(ii)),L);
    fprintf('Ok (%gs)\n',toc);
end
tFDM = toc(tstart);

% Do each r value
Tr = [];
Hr = [];
tic;
for ii = 1:length(r)
    fprintf('Split %u/%u: r=%g... ',ii,length(r),r(ii)); 
     G = zeros(N,INTERP_COLS);
     for ij = 1:INTERP_COLS
         G(:,ij) = (A{ij}(ii,:));
     end
    
    [Tr1 Hr1 FLAG RELRES ITER] = splitTPH(G,ip);
    Tr = [Tr Tr1];
    Hr = [Hr Hr1];
    fprintf('%u iters (res=%g)\n',ITER,RELRES);
end

disp('====================');
fprintf('Total FDM time=%gs\n',tFDM);
fprintf('Total split time=%gs\n',toc);

% Prepare the guide vectors.
r(1) = 0; % smudge the first r point to make it zero
zh = linspace(2*zp(1),2*zp(N),2*N-1);
zt = linspace(0,zp(N)-zp(1),N);

end

function [t h FLAG RELRES ITER] = splitTPH(Gcols,ip)
%SPLITTPH Splits the TPH matrix G into T and H
%   Symmetry along the main diagonal is assumed (reciprocity statement)
%   Gcols is the column of G(z,zp(ip)) evaluated at each zp value

[M n_ip] = size(Gcols);
assert(n_ip==numel(ip), 'Input must be ordered such that Gcols has as many columns as ip has elements.')

t=zeros(2*M-1,1);
h=zeros(2*M-1,1);

% Construct the linear algebra matrix
A = [];
for ii = 1:n_ip
    A = [A; get_tn_sym(ip(ii)), get_hn(ip(ii))];
end


% Construct the b vector
b = reshape(Gcols,n_ip*M,1);

% Solve using LSQR
[x,FLAG,RELRES,ITER] = lsqr(A,b,1e-12,500);

% Extract
t = x(1:M);
h = x(M+1:end);

    function Hn = get_hn(ii)
    %GETHN Gets the M x (2M-1) coefficient matrix for the coefficients of H
    %  for a particular index
        Hn = [sparse(M,ii-1), speye(M), sparse(M,M-ii)];
    end

    function Tn = get_tn_sym(ii)
    %GETHN Gets the M x M coefficient matrix for the symmetric T of 
    % a particular index
        %Tn = [sparse(M-ii,ii-1), speye(M-ii), sparse(M-ii,1);
        %    sparse(ii,M-ii),spjay(ii)];
        Tn = [spjay(ii),sparse(ii,M-ii);
            sparse(M-ii,1), speye(M-ii), sparse(M-ii,ii-1)];
    end
     
    function out = spjay(N)
    %SPJAY like speye, except for the antidiagonal identity J.
        out = sparse(1:N,N:-1:1,ones(1,N));
    end
end
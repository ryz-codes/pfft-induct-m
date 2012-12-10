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
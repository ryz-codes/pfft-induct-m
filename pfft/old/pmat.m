function [P] = pmat(s,qpnt,qwei,qexc)
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
%

persistent invG

ord = s.L; % Order is equal to the length of the stencil
sten_coord = s.getStenCoord; % retrieve coords
Ns = size(sten_coord,1); % number of stencil points

dims = size(qpnt,2);

qpnt = cell2mat(qpnt(:));
Nq = size(weights,2);
Nb = size(qpnt,1)/Nq; % number of bases 

% Get inversion of G
if ~any(any(invG))
    invG = pinv(f_r(ord,sten_coord));
end

% initialize

% Iterate and quadrature each basis function
do_wei = (nargin > 2) && any(qwei);
do_exc = (nargin > 3) && any(qexc);
% for k = 1:Nb
%     this_qpnt = qpnt(k,:);
%     %this_qpnt = reshape(this_qpnt,Nq,dims);
%     this_line = f_r(ord,this_qpnt)*invG;
%     
%     if do_wei
%         this_line = qwei*this_line;
%     end
%     if do_exc
%         this_line = qexc(k)*this_line;
%     end
%     P(:,k) = this_line;
% end

P = (f_r(ord,qpnt)*invG).';

if do_exc
    bsxfun(@times,P,qexc.')
end
end

function [out] = f_r(ord, pnts)
%F_R Polynomial interpolation function
%   Gives a row vector of the interpolation basis functions
%   pnts must be column vectors, one column for each dimension.
persistent s
persistent prev_ord

[num_ele dims] = size(pnts); 
if ~any(any(s)) || (ord ~= prev_ord) || (size(s,1)~=dims) % we don't have s
    % make s
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
    s = [];
    for ii = 1:dims
        temp = s_all(ii,:);
        s(ii,:) = temp(s_sel);
    end
    
    % Store
    prev_ord = ord;
end

% Multiply and generate f_r
num_exp = size(s,2);
out = 1;
for ii = 1:dims
    out = out .* (repmat(pnts(:,ii),1,num_exp))...
        .^(repmat(s(ii,:),num_ele,1));
end
end



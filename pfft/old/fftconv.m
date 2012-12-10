function [ field ] = fftconv( fftkernel, charge)
%FFTCONV Convolves in all possible dimensions
%   fftkernel should be generated by g.getCoordBox and have dimensions
%   equal to the dimension of the convolution. Charge should have
%   dimensions exactly half that of fftkernel. Charge can either be a
%   vector or a matrix block

% Size of the box / cube considered
ii = 1;
this_dim = size(fftkernel,1);
while this_dim > 1
    siz(1,ii) = this_dim/2;
    ii = ii +1;
    this_dim = size(fftkernel,ii);
end

num_ele = prod(siz);
assert(numel(charge)==num_ele,'Number of elements between kernel and charge do not match!');

vec = false; 
if size(charge,1) == num_ele && ...
        numel(siz)>1% was charge inputted as a vector?
    vec = true;
    charge = reshape(charge,siz);
end

% Pad zeros
[charge idx] = ConstantPad(charge,siz);
% Convolve
this_proj = ifftn(fftn(charge).*fftkernel);
% Remove the extra rows
field = real(this_proj(idx{:}));

if vec
    field = real(field(:));
end
end

%%%
%%% ConstantPad
%%%
function [b idx] = ConstantPad(a, padSize)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
sizeB = zeros(1,numDims);
for k = 1:numDims
    M = size(a,k);
    idx{k}   = 1:M;
    sizeB(k) = M + padSize(k);

end

% Initialize output array with the padding value. 
if numel(sizeB) == 1
    b         = zeros(sizeB,1);
else
b         = zeros(sizeB);
end
b(idx{:}) = a;
end
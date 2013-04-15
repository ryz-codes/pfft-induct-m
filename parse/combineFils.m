function [ fils ] = combineFils( varargin)
%COMBINEFILS Summary of this function goes here
%   Detailed explanation goes here
fils = cell(1,4);
for ii = 1:nargin
    temp = varargin{ii};
    assert(numel(temp)==4);
    for ij = 1:4
        fils{ij} = [fils{ij};temp{ij}];
    end
end

end


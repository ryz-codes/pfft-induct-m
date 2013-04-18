function [ g_cell ] = sweep_tph(w,id)
%SWEEP_TPH Summary of this function goes here
%   This function only sweeps defaultL(1) and defaultL(5)
if nargin == 1
    id = 1;
end
if id == 1 
    bnds = {0.7,0.999};
elseif id == 5
    bnds = {0.01,0.99};
else
    error('Only id == 1 or 5 supported!!');
end

N = length(w);
L = defaultL(id);

bnds3_scale = L.bnds(3) * sqrt(L.w); % scaling factor for half-space
g_cell=cell(N,1);
for ii = 1:N
    L.w = w(ii);
    
    if id == 1 %Only in the half-space case
        % adjust thickness of the halfspace layer to be 5 skin depths
        L.bnds(3) = bnds3_scale / sqrt(w(ii)); 
    end
    
    g_cell{ii} = tph(L,bnds{:});
end


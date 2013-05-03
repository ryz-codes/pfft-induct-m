function [O,L,W,H] = path2segs(path,rad)
%SEGSFROMPATH Summary of this function goes here
%   Detailed explanation goes here
pxrad = 4;
N = size(path,1)-1;
O = cell(N,1); L = cell(N,1); W = cell(N,1); H = cell(N,1);
for ii = 1:N
    [O{ii},L{ii},W{ii},H{ii}] = seg2fils(path(ii,:),path(ii+1,:),rad,pxrad);
end

end


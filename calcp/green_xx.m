function [G] = green_xx(r,z,zp)
%GREEN Green's function wrapper for Gxx

persistent g

if isempty(g)
    %addpath([pwd '/../green']) 
    load default4.mat;
end

G = g.lookup(r,z,zp);
end


%TPH_WRAP
% A wrapper of the tph green's function class for the pfft class solver.
% in - input of coordinates, with in{1} = x, in{2} = y, in{3} = z-zp or
%      z+zp
% type - either 'T' or 'H'
function out = tph_wrap(g,type,in)
    r = sqrt(in{1}.^2 + in{2}.^2);
    z = in{3};
    if strcmp(type,'T')
        out = g.lookupT(r,z,0);
    elseif strcmp(type,'H');
        out = g.lookupH(r,z,0);
    else
        error('type is either T or H');
    end
end
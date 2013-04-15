function plotfil(from, to)
%PLOTFIL Summary of this function goes here
%   Detailed explanation goes here
if(size(from,1) == 3);

    % strip z axis
    from = from(1:2,:);
    to = to(1:2,:);
elseif size(from,1) ~= 2
    error('fail!!');
end

% MATLAB plot function plots 
fromx = from(1,:);
fromy = from(2,:);
tox = to(1,:);
toy = to(2,:);

plot([fromx;tox],[fromy;toy]);


end


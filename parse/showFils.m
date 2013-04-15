function [ output_args ] = showFils(O,L,W,H)
%SHOWFILS Summary of this function goes here
%   Detailed explanation goes here
%       7----8
%      /|   /|
%     / 5----6
%    / /  / /
%   / /  / /
%  3-/--4 /
%  |/   |/
%  1----2
%
% Make eight vertices
%
% If the number of fils to draw exceeds 200, we won't draw the front and 
% back faces.

draw_face = size(O,1) < 200;
if nargin == 1 && numel(O) == 4 %if input is a cell
    L = O{2};
    W = O{3};
    H = O{4};
    O = O{1};
end

vert{1} = O;
vert{2} = O+W;
vert{3} = O+H;
vert{4} = O+W+H;
vert{5} = O+L;
vert{6} = O+L+W;
vert{7} = O+L+H;
vert{8} = O+L+W+H;

% Make six faces
if draw_face
    [f1x f1y f1z] = getface (1,2,4,3);
    [f6x f6y f6z] = getface (5,6,8,7);
end
[f2x f2y f2z] = getface (1,5,7,3);
[f3x f3y f3z] = getface (2,6,8,4);
[f4x f4y f4z] = getface (3,4,8,7);
[f5x f5y f5z] = getface (1,2,6,5);


% Draw the faces
clf;
for ii = 1:size(O,1)
    if draw_face
        patch(f1x(ii,:), f1y(ii,:), f1z(ii,:), [1 0.5 0.5]); %pink
        patch(f6x(ii,:), f6y(ii,:), f6z(ii,:), [1 0.5 0.5]); %pink
    end
    patch(f2x(ii,:), f2y(ii,:), f2z(ii,:), 'red');
    patch(f3x(ii,:), f3y(ii,:), f3z(ii,:), 'red');
    patch(f4x(ii,:), f4y(ii,:), f4z(ii,:), 'red');
    patch(f5x(ii,:), f5y(ii,:), f5z(ii,:), 'red');
end
view(3);
axis('equal');

    % Private function to sort out the vertices
    function [fx, fy, fz] = getface (a,b,c,d)
        fx = [vert{a}(:,1) vert{b}(:,1) vert{c}(:,1) vert{d}(:,1)];
        fy = [vert{a}(:,2) vert{b}(:,2) vert{c}(:,2) vert{d}(:,2)];
        fz = [vert{a}(:,3) vert{b}(:,3) vert{c}(:,3) vert{d}(:,3)];
    end
end

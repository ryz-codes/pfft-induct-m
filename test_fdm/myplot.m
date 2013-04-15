function myplot(h,r1,G1,r2,G2,offset)
%MYPLOT Summary of this function goes here
%   Detailed explanation goes here
if nargin == 5, offset=0; end

num_items = size(G1,2);

c_set =      [0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500
         0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500];
m_set = {'o','p','v','^','s','.','+','d','*','>','<','x','h'};

hold(h,'on');


for ii = 1:num_items
    % Plot markers
    plot(h,r1,G1(:,ii),'color',c_set(ii+offset,:),'marker',m_set{ii+offset},...
        'linestyle','none','MarkerSize',6);
end

% Plot lines
for ii = 1:num_items
    plot(h,r2,G2(:,ii),'color',c_set(ii+offset,:));
end

% Hold off
hold(h,'off');

end


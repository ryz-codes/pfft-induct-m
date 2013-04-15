function M = filquad(O,L,W,H)

[nodes, weights] = nwspgr('KPU', 4, 6);

M = 0;
% Make a quadrature in each cell.
for ii = 1:numel(weights)
    offset1 = W(:,1:3)*nodes(ii,1) + ...
        H(:,1:3)*nodes(ii,2);
    offset2 = W(:,4:6)*nodes(ii,3) + ...
        H(:,4:6)*nodes(ii,4);
    
    from1 = O(:,1:3) + offset1;
    to1 = O(:,1:3) + L(:,1:3) + offset1;
    
    from2 = O(:,4:6) + offset2;
    to2 = O(:,4:6) + L(:,4:6) + offset2;
    
    M = M + mutualfil(from1',to1',from2',to2') * weights(ii);
end
disp(numel(weights))
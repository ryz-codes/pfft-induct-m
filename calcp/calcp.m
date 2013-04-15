function [V M] = calcp(I,O,L,W,H,g_xx)
%CALCP Calculates the voltage across filaments caused by the Gxx subtracted
%Green's function
% This method is vectorized, and works best when used in a vectorized
% fashion.
%
% Inputs for m filaments self interactions:
%    I - current vector, m x 1 column vector
%    O - origin vector, m x 2 or m x 3 column vector
%    L - length vector, m x 2 or m x 3 column vector
%    W - width vector, m x 3 column vector
%    H - height vector, m x 3 column vector

    kronpairs = @(a,N) [kron(a,ones(N,1)) kron(ones(N,1),a)];

    % Collect information and throw exceptions if input is problematic
    [N shouldbeone] = size(I); % number of filaments
    assert(shouldbeone == 1,'Current vector (I) should only have one column!');
    
    % To do:
    %   - 2D filaments
    %   - zero volume filaments
    %   - source to field interactions
    
    % Generate interaction pairs. This part is obviously O(N^2), shown by
    % the use of the kronicker product.
    Op = kronpairs(O,N);
    Lp = kronpairs(L,N);
    Wp = kronpairs(W,N);
    Hp = kronpairs(H,N);

    % Matlab quadrature the subtracted component
    if nargin == 6
        Mxx = mf6_xx(Op, Lp, Wp, Hp, g_xx);
    else
        Mxx = mf6_xx(Op, Lp, Wp, Hp);
    end
    Mxx = reshape(Mxx,N,N);
    
    % Freespace component is done by FastHenry port
    M = mutual(Op.', Lp.', Wp.', Hp.'); 
    M = reshape(M,N,N);

    V=(M+Mxx)*I;
    %V=Mxx*I;
end
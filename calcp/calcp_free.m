function [V M] = calcp_free(I,O,L,W,H)
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

    % Collect information and throw exceptions if input is problematic
    [N shouldbeone] = size(I); % number of filaments
    assert(shouldbeone == 1,'Current vector (I) should only have one column!');
    
    % To do:
    %   - 2D filaments
    %   - zero volume filaments
    %   - source to field interactions
    
    % Generate interaction pairs. This part is obviously O(N^2), shown by
    % the use of the kronicker product.
    if size(O,2) == 3
        Op = kronpairs(O,N);
        Lp = kronpairs(L,N);
        Wp = kronpairs(W,N);
        Hp = kronpairs(H,N);
    elseif size(O,2) == 6
        Op = O;
        Lp = L;
        Wp = W;
        Hp = H;
    else
        error('Only three dimensional systems implemented');
    end
        

    % Matlab quadrature the subtracted component
    
    % Freespace component is done by FastHenry port
    M = mutual(Op.', Lp.', Wp.', Hp.'); 
    M = reshape(M,N,N);

    V=(M)*I;
end


function c = kronpairs(a,N)
    c=[kron(a,ones(N,1)) kron(ones(N,1),a)];
end

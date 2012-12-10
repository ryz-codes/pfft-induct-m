function test_integrand
ang = pi;
bnd_a = -1;
bnd_b = 1;

%a = linspace(-10,10,1000);
a = logspace(-10,1);
M = gen_M(a,ang);
figure(2);semilogx(a,M);

% Direct quadrature
[result ebnd] = quadgk(@(x) gen_M(x,ang), -1, 1, 'RelTol', 1e-10)



end

function M = gen_M(a,ang)

siz = size(a);
a = reshape(a,1,[]);
a = padarray(a,[1,0]);

N = length(a);

from1 = [0 0 0]';
to1 = [1 0 0]';

% Get
from2 = bsxfun(@plus,from1,a);
to2 = [cos(ang) sin(ang) 0]';
to2 = bsxfun(@plus,to2,a);
from1 = repmat(from1,1,N);
to1 = repmat(to1,1,N);

M = mutualfil(from1,to1,from2,to2);
M = reshape(M,siz);
end
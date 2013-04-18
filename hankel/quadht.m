function [c, errest] = quadht(fun, a, k, nu, quadord)
%QUADHT Evaluate a hankel transform integral using Clenshaw-Curtis
%Quadrature.
%   The method included here evaluates the hankel integral from 0 to a:
%      [integral from 0 to a] fun(k) * besselj(nu,k*r) * r * dr 
%   Using the Clenshaw-Curtis quadrature. 
%
%   The quadrature replaces the following code
%      [A,r,k] = fht(fun,R,K,nu);
%   with
%      A = quadht(fun,R,k,nu);
%
%   or in reverse:
%      [~,r,k,I] = fht([],R,K,nu);
%   with
%      A = quadht(fun,K,r,nu);
persistent x
persistent w
persistent w2
if nargin < 5, quadord=1e5;end

% Calculate quadrature nodes and weights if they haven't been calculated
if mod(quadord,2) == 0 % make odd
    quadord = quadord +1;
end
quad2ord = (quadord+1)/2;
[x,w] = fclencurt(quadord,0,1);
[~,w2] = fclencurt(quad2ord,0,1);

siz = size(k);
assert(length(k)==numel(k)); % must be a vector 

% Perform quadrature
r = x*a;
k = k(:).'; % row vector
Jkr = besselj(nu,bsxfun(@times,k,r)); % matrix
rf = r.*fun(r);
integrand = bsxfun(@times,Jkr,rf);
integrand(isnan(integrand)) = 0; %set all nan to zero.
c = a*w.'*integrand;

% Quadrature lower order.
integrand = integrand(1:2:end,:);
c2 = a*w2.'*integrand;

% Error is infinity (aka max) norm of the differences.
errest = norm(c-c2,inf);

c = reshape(c,siz);
end


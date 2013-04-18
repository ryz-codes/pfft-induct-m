function [c] = endcor(fun, a, k, nu, quadord)
%ENDCOR QFHT end correction term via Gaussian Quadrature 
%   A perfect hankel integral goes from 0 to infinity. Unfortunately, in
%   practice this integral is only evaluated from a to b, where a
%   approximates 0 and b approximates infinity.
%
%   The method included here evaluates the hankel integral from 0 to a:
%      [integral from 0 to a] fun(k) * besselj(nu,k*r) * r * dr 
%   and thus gives an estimate for trunctation error at a. With this, the
%   error is governed only by the upper bound b.
%
%   The following code can be combined with the fht or ifht code:
%      [A,r,k] = fht(fun,R,K,nu);
%      A = A + endcor(fun,r(1),k,nu);
%
%   or in reverse:
%      [~,r,k,I] = fht([],R,K,nu);
%      A = ifht(fun(k),k,r,I) + endcor(fun,k(1),r,nu);
%
%   LIMITATIONS: the function fun must be smooth between 0 and a (<20th 
%   order polynomial).
%
%   Inspired by not directly equivalent to the method described by:
%       G. Agrawal and M. Lax, "End correction in the quasi-fast Hankel 
%       transform for optical propagation problems," Opt. Lett.  6, 
%       171-173 (1981).

if nargin < 5, quadord=20;end

% Gaussian quadrature nodes and weights for order = 20
[x,w] = fclencurt(quadord,0,1);
% x =[0.003435700407453
%    0.018014036361043
%    0.043882785874337
%    0.080441514088891
%    0.126834046769925
%    0.181973159636742
%    0.244566499024586
%    0.313146955642290
%    0.386107074429177
%    0.461736739433251
%    0.538263260566749
%    0.613892925570823
%    0.686853044357710
%    0.755433500975414
%    0.818026840363258
%    0.873165953230075
%    0.919558485911109
%    0.956117214125663
%    0.981985963638957
%    0.996564299592547];
% 
% 
% w =[0.008807003569575
%    0.020300714900194
%    0.031336024167055
%    0.041638370788352
%    0.050965059908620
%    0.059097265980759
%    0.065844319224588
%    0.071048054659191
%    0.074586493236302
%    0.076376693565363
%    0.076376693565363
%    0.074586493236302
%    0.071048054659191
%    0.065844319224588
%    0.059097265980759
%    0.050965059908620
%    0.041638370788352
%    0.031336024167055
%    0.020300714900194
%    0.008807003569575];

siz = size(k);
assert(length(k)==numel(k)); % must be a vector 

r = x*a;
k = k(:).'; % row vector
Jkr = besselj(nu,bsxfun(@times,k,r)); % matrix
rf = r.*fun(r);
integrand = bsxfun(@times,Jkr,rf);
integrand(isnan(integrand)) = 0; %set all nan to zero.
c = a*w.'*integrand;
c = reshape(c,siz);
end


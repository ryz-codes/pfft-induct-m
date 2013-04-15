function [L] = analy_coils_free(R,c)
% Calculates the free-space inductance. Assumes the spiral turns can be
% approximated as circles. Uses the hankel integral for mutual terms and
% the rectangular form for self.
% 
persistent k1
persistent r1
persistent I
%% Prepare Kernel
if isempty(I)
    [~,k1,r1,I]=fht(@(x) 1,max(R)*10,2*pi*10/min(R),1);
end

%% Preparation
addpath([pwd '/hankel']) % Include the package

mu=1e-7*4*pi;

% Setup analytical Green's function
Ar_ker = mu * pi ./ k1; 

% Superposition the contribution of the rings
L = 0;
for ii = 1:length(R)
    % self term
    L = L + analySqRing(R(ii),c);
    
    % mutual term
    R1 = R(R~=R(ii));
    Ar_ker_num = besselj(1,bsxfun(@times,k1(:).',R1(:)));
    Ar_ker_num = Ar_ker.*sum(bsxfun(@times,Ar_ker_num,R1(:)),1);
    A = r1.*ifht(Ar_ker_num,k1,r1,I);
    L = L + interp1(r1,A,R(ii));
    
    ii
end

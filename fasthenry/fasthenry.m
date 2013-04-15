function [Zt] = fasthenry(path,rad)
%FASTHENRY Terminal impedance extraction akin to FastHenry. Only difference
%is that only a path is given, and all conductors are assumed to be round.

% path is either a matrix or a cell array. 
% rad is the radius of the path. 
path = linspace(0,2*pi).';
path = [cos(path),sin(path),zeros(size(path,1),1)];
rad = 1e-3;
w = 2*pi*logspace(1,8,15);

[O,L,W,H] = path2segs(path,rad);

% Make segment-based mesh matrix
Nfils = size(O{1},1); % number of fils per segment
Nsegs = size(O,1);
M1 = sparse(1,(0:Nsegs-1)*Nfils+1,1,1,Nfils*Nsegs);
M2 = [ones(Nfils-1,1), spdiags(-ones(Nfils-1,1),0,Nfils-1,Nfils-1)];
M3 = cell(1,Nfils);
for ii = 1:Nsegs
    M3{ii} = M2;
end
M = [M1; blkdiag(M3{:})];
Nmesh = size(M,1);

% Expand out the branches
O = vertcat(O{:}); L = vertcat(L{:}); W = vertcat(W{:}); H = vertcat(H{:});

% Set up the grid
induct(O,O+L+W+H);
induct(ones(Nfils*Nsegs,1),O,L,W,H);

% Objective vector
Vm0 = zeros(Nmesh,1); Vm0(1) = 1;

% Solve for each frequency
for ii = 1:length(w)
    Im0 = gmres(@MZMT,Vm0,[],1e-3,30);

    Zt(ii)=1/Im0(1);
end

% Plot cross-sectional current density
figure(1);
subplot(121);semilogx(w/2*pi,imag(Zt)./w); title('Inductance vs f');
subplot(122);loglog(w/2*pi,real(Zt)); title('Resistance vs f');

    % Multiply by mesh impedance matrix
    function Vm = MZMT(Im)
        Ib = M'*Im;
        Vm = M*(1j*w(ii)*induct(Ib)+resist(Ib,[],L,W,H));
    end
end
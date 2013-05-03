function bwmap
%BWMAP Blue, red white colormap suitable for black and white and color
%printing
N = 10;
N2 = N/2;

a = linspace(0,1,N)';
c = [zeros(N/2,2),a(N2+1:end); % deep blue to blue
     a*[1 1],ones(N,1)]; % blue to white
c = [c;rot90(c,2)];
colormap(c)
 

end


%%% Tests the fht files for some common transform pairs.

%% Test the forward-backward transform error
spfun = @(r,v) (r.^v).*(r<=1);
fqfun = @(k,v) besselj(v+1,k)./k;

ORD=1;
[F,k,r,I] = fht(@(x)spfun(x,ORD),2,100,ORD);
f = ifht(F,k,r,I);
subplot(211)
plot(r,f,r,spfun(r,ORD));
title('Spatial');
subplot(212)
semilogx(k,F,k,fqfun(k,ORD));
title('Frequency');

%% Test the backward transform
ORD = 0;
R = 2;
K = 100;
[~,k1,r1,I1]=fht([],R,K,ORD,4,6);

% Setup reverse transform function
subplot(311)
loglog(k1,abs(fqfun(k1,ORD)));

% Perform reverse transform
f1 = ifht(fqfun(k1,ORD),k1,r1,I1)+ ...
    endcor(@(x)fqfun(x,ORD),k1(1),r1,ORD); 
r2 = linspace(r1(1),r1(end),10); % pick 10 points to test.
f2 = endcor(@(x)fqfun(x,ORD),k1(end)*1e2,r2,ORD,5e5); % use an intense quadrature to evaluate
f3 = quadiht(@(x)fqfun(x,ORD),k1(end)*1e2,r2,ORD);


% Plot compare to actual values
f0 = spfun(r1,ORD); % actual values
subplot(312)
plot(r1,f0,r1,f1,r2,f2);

subplot(313)
er = median(abs(f1(r1<0.5)-f0(r1<0.5)));
plot(r1(r1<0.5),f1(r1<0.5)-f0(r1<0.5))

%% Error estimates
A0 = trapz(r1,f0);
A1 = trapz(r1,f1);

err1 = (A1-A0)/A0 % actual error
e_est = max(abs(interp1(r1,f1,r2)-f2))



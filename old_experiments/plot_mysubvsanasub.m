addpath([pwd '/test_fdm'])

load new_pfft;
% Get the L and R
L_c = real(Vsub);
L_a = real(Vsuba);
L_e = abs(L_c-L_a);
R_c = -imag(Vsub).*w(:);
R_a = -imag(Vsuba).*w(:);

figure(1)
myplot(gca,w,L_c,w,L_a);
set(gca,'XScale','log')
xlim('auto')
ylabel('\DeltaL [H]');
xlabel('w [rads^{-1}]');
legend('This work','Analy. sol.');
grid on

figure(2)
myplot(gca,w,R_c,w,R_a);
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim('auto')
ylabel('\DeltaR [H]');
xlabel('w [rads^{-1}]');
legend('This work','Analy. sol.');
grid on
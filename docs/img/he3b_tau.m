#!/usr/bin/octave-cli -f
pkg load he3lib

% Plot quasiparticle lifetime

figure;
hold on;

p1=0;
p2=30;

ttc=0.15:0.002:1.4;

subplot(1,2,1); hold on;
semilogy(ttc, he3_tau_n0(ttc,p1),   'r-');
semilogy(ttc, he3_tau0(ttc, p1),    'r-', 'linewidth', 2);
semilogy(ttc, he3_tau_n_av(ttc,p1), 'b-');
semilogy(ttc, he3_tau_av(ttc, p1),  'b-', 'linewidth', 2);
semilogy(ttc, he3_tau0lt(ttc,p1),   'g-');

semilogy(ttc, he3_tau_n0(ttc,p2),   'r--');
semilogy(ttc, he3_tau0(ttc, p2),    'r--', 'linewidth', 2);
semilogy(ttc, he3_tau_n_av(ttc,p2), 'b--');
semilogy(ttc, he3_tau_av(ttc, p2),  'b--', 'linewidth', 2);
semilogy(ttc, he3_tau0lt(ttc,p2),   'g--');

xlabel('T/Tc')
ylabel('\tau, s')
legend('he3\_tau\_n0', 'he3\_tau0', 'he3\_tau\_n\_av', 'he3\_tau\_av', 'he3\_tau0lt')
title('quasiparticle lifetime (0bar, 30bar)')
text(1,2e-6, '0bar')
text(1,1e-7, '30bar')
ttc=0.15:0.002:1.4;

subplot(1,2,2); hold on;
l1=he3_vf(p1)*he3_tau_n0(ttc, p1);
l2=he3_vf(p2)*he3_tau_n0(ttc, p2);
semilogy(ttc, l1, 'r-');
semilogy(ttc, he3_fpath(ttc,p1),                'r-', 'linewidth', 2);
semilogy(ttc, l2, 'r--');
semilogy(ttc, he3_fpath(ttc,p2),            'r--', 'linewidth', 2);

text(1,1e-2, '0bar')
text(1,3e-4, '30bar')

xlabel('T/Tc')
ylabel('l, cm')
legend('he3\_vf*he3\_tau\_n0', 'he3\_fpath')
title('quasiparticle mean free path (0bar, 30bar)')


print -dpng he3b_tau.png "-S800,300" "-F:6"

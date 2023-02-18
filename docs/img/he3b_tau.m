#!/usr/bin/octave-cli -f
pkg load he3lib

% Plot quasiparticle lifetime

figure;
hold on;

p=0;

ttc=0.15:0.01:1.4;

subplot(1,2,1); hold on;
semilogy(ttc, he3_tau_n0(ttc,p),   'b--');
semilogy(ttc, he3_tau_n_av(ttc,p), 'r--');
semilogy(ttc, he3_tau0(ttc, p),    'b-');
semilogy(ttc, he3_tau_av(ttc, p),  'r-');

xlabel('T/Tc')
ylabel('\tau, s')
legend('he3\_tau\_n0', 'he3\_tau\_n\_av', 'he3\_tau0', 'he3\_tau\_av')
title('quasiparticle lifitime (P=0bar)')

ttc=0.15:0.01:1.4;
tn0=he3_tau_n0(ttc, p);
vf=he3_vf(p);

subplot(1,2,2); hold on;
semilogy(ttc, he3_fpath(ttc,p),      'r-');
semilogy(ttc, he3_visc_fpath(ttc,p), 'b-');
semilogy(ttc, vf*tn0,               'r--');
xlabel('T/Tc')
ylabel('l, cm')
legend('he3\_fpath', 'he3\_visc\_fpath', 'he3\_vf()*he3\_tau\_n0()')
title('quasiparticle mean free path (P=0bar)')


print -dpng he3b_tau.png "-S800,300"

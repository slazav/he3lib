#!/usr/bin/octave-cli -f
pkg load he3lib

figure; clf; hold on;
ttc = 0.1:0.002:1.2;
p=21;

tn0=he3_tau_n0(ttc, p);
vf=he3_vf(p);
gap=he3_trivgap(ttc,p);

l0=tn0.*vf;

% table 2 in Einzel-1987
ttc1 = [1.0 0.6 0.5 0.4 0.35 0.3 0.25 0.2 0.15 0.1];
mfp1 = [1.70e-4 4.68e-4 9.13e-4 2.42e-3 4.78e-3 1.17e-2 3.99e-2 2.45e-1 4.87 1.82e3];

%  semilogy(ttc, he3_tau_av(ttc, p) , 'r.-');
%  semilogy(ttc, he3_tau_n_av(ttc, p) , 'b.-');

semilogy(ttc, he3_fpath(ttc, p) , 'r.-');
semilogy(ttc, he3_fpath_n(ttc, p) , 'b.-');
semilogy(ttc, tn0.*vf , 'm.-');
semilogy(ttc1, mfp1 , 'g*');

print -dpng -color plot_mfp.png

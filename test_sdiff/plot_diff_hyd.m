#!/usr/bin/octave -qf
  % Compare \omega=0 case with hydrodinamic formula
  addpath ../matlab

  figure; clf; hold on;
  ttc = [0.2:0.01:0.95 0.951:0.001:1 1.01:0.01 1.5];
  p=30;

  f=5e5;

plot(ttc, he3_diff_par_zz(ttc, p, 0), 'r-', 'linewidth',2);
plot(ttc, he3_diff_hpar_zz(ttc, p),   'b-');

plot(ttc, he3_diff_perp_zz(ttc, p, 0), 'm-', 'linewidth',2);
plot(ttc, he3_diff_hperp_zz(ttc, p),   'c-');

  xlim([0.2 1.5])
  ylim([0 0.1])
  xlabel('T/T_c')
  ylabel('D, cm^2/s')

  print plot_diff_hyd.eps -deps -color "-S800,600"

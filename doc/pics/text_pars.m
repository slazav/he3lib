#!/usr/bin/octave -qf
  addpath ../../matlab

  figure; clf; hold on;

  subplot(3,2,1); hold on;
  plot_tpdep(@he3_text_a, 'northwest');
  title('a, erg/cm^3 1/G^2 vs T/T_c');

  subplot(3,2,2); hold on;
  plot_tpdep(@he3_text_ldv);
  title('\lambda_{DV}, erg/cm^3 1/(cm/s)^2 vs T/T_c');
  ylim([0 2.5e-7])

  subplot(3,2,3); hold on;
  plot_tpdep(@he3_text_lhv);
  title('\lambda_{HV}, erg/cm^3 1/(G cm/s)^2 vs T/T_c');

  subplot(3,2,4); hold on;
  plot_tpdep(@he3_text_d);
  title('d, erg/cm^2 1/G^2 vs T/T_c');

  subplot(3,2,5); hold on;
  f = @(ttc,p) (he3_text_llh(ttc,p,1));
  plot_tpdep(f);
  title('\lambda_{LH} at \Omega=1 rad/s vs T/T_c');

  subplot(3,2,6); hold on;
  f = @(ttc,p) he3_text_lo(ttc,p,1);
  plot_tpdep(f);
  title('\lambda/\Omega vs T/T_c');

  print text_pars.eps -deps -color "-S1200,1600"

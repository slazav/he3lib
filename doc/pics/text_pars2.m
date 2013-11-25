#!/usr/bin/octave -qf
  addpath ../../matlab

  figure; clf; hold on;

  subplot(2,2,1); hold on;
  plot_tpdep(@he3_text_d, 'southwest');
  title('d, erg/cm^2 1/G^2 vs T/T_c');

  f = @(ttc,p) (he3_text_llh(ttc,p,1));
  subplot(2,2,2); hold on;
  plot_tpdep(f);
  title('\lambda_{LH} at \Omega=1 rad/s vs T/T_c');

  f = @(ttc,p) (5/2.0 * he3_text_llh(ttc,p,1)./he3_text_a(ttc,p));
  subplot(2,2,3); hold on;
  plot_tpdep(f);
  title('\lambda/\Omega vs T/T_c');

  print text_pars2.eps -deps -color "-S1200,800"

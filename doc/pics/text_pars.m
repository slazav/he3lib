#!/usr/bin/octave -qf
  addpath ../../matlab

  figure; clf; hold on;

  subplot(4,2,1); hold on;
  plot_tpdep(@he3_text_a, 'northwest');
  title('a, erg/cm^3 1/G^2 vs T/T_c');

  subplot(4,2,2); hold on;
  plot_tpdep(@he3_text_ldv);
  title('\lambda_{DV}, erg/cm^3 1/(cm/s)^2 vs T/T_c');
  ylim([0 2.5e-7])

  subplot(4,2,3); hold on;
  plot_tpdep(@he3_text_lhv);
  title('\lambda_{HV}, erg/cm^3 1/(G cm/s)^2 vs T/T_c');

  subplot(4,2,4); hold on;
  plot_tpdep(@he3_text_lg2);
  title('\lambda_{G2}, erg/cm vs T/T_c');
  ylim([0 3e-11])

  subplot(4,2,5); hold on;
  plot_tpdep(@he3_text_delta);
  title('\delta vs T/T_c');
  ylim([-0.35 0])

  subplot(4,2,6); hold on;
  plot_tpdep(@he3_text_vd);
  title('v_d, cm/s vs T/T_c');
  ylim([0 0.08])

  subplot(4,2,7); hold on;
  plot_tpdep(@he3_text_xid);
  title('\xi_d, cm vs T/T_c');

  subplot(4,2,8); hold on;
  f=@(t,p) (he3_text_xih(t,p, 250));
  plot_tpdep(f);
  title('\xi_h at 250 Oe, cm vs T/T_c');

  print text_pars.eps -deps -color "-S1200,1600"

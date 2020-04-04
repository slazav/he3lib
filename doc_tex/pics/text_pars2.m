#!/usr/bin/octave -qf
  graphics_toolkit("gnuplot")
  figure; clf; hold on;


  subplot(2,2,1); hold on;
  plot_tpdep(@he3_text_xid);
  title('\xi_d, cm vs T/T_c');

  subplot(2,2,2); hold on;
  f=@(t,p) (he3_text_xih(t,p, 250));
  plot_tpdep(f);
  title('\xi_h at 250 Oe, cm vs T/T_c');

  subplot(2,2,3); hold on;
  plot_tpdep(@he3_text_vd);
  title('v_d, cm/s vs T/T_c');
  ylim([0 0.08])

  print text_pars2.eps -deps -color "-S500,400"

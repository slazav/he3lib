#!/usr/bin/octave-cli -qf

  pkg load he3lib
  figure; clf; hold on;

  subplot(1,2,1); hold on;
  plot_tpdep(@he3_cperp);
%  plot_tpdep(@he3_cpar);
  title('c_{perp}, cm/c vs T/T_c');
%  f=@(ttc,p)(1100*sqrt(1.5)./he3_vf(34.3).*he3_vf(p).*sqrt(1-ttc) );
%  plot_tpdep(f);
  ylim([0 2600])

  subplot(1,2,2); hold on;
  f = @(ttc,p) (he3_cpar(ttc,p));
  plot_tpdep(f);
  title('c_{parallel}, cm/s vs T/T_c');
  ylim([0 2600])

  print text_gr.eps -deps "-S500,200"

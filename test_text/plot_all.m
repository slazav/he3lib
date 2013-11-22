#!/usr/bin/octave -qf
  addpath ~/he3lib/lib/matlab
  plot_tpdep(@he3_text_a);
  title('textural parameter a');
  ylabel('a');

  plot_tpdep(@he3_text_delta);
  title('textural parameter \delta');
  ylabel('\delta');

  plot_tpdep(@he3_text_lhv);
  title('textural parameter \lambda_{HV}');
  ylabel('\lambda_{HV}');

  plot_tpdep(@he3_text_lg2);
  title('textural parameter \lambda_{G2}');
  ylabel('\lambda_{G2}');

  plot_tpdep(@he3_text_vd);
  title('textural parameter v_d');
  ylabel('v_d');
  ylim([0 200])

  plot_tpdep(@he3_text_ksid);
  title('textural parameter \ksi_d');
  ylabel('\ksi_d');

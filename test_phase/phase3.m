#!/usr/bin/octave -qf

  addpath ../matlab

  figure; clf; hold on;
  t=0.5:0.01:3;
  plot(t, he3_pmelt(t/1000), 'b-');

  plot(t, he3_pmelt_plts(t/1000), 'r-');
  plot(he3_ta_plts, he3_pa_plts, 'ro');
  plot(he3_tb_plts, he3_pb_plts, 'ro');
  plot(he3_ts_plts, he3_ps_plts, 'ro');


  p=34.32:0.001:34.4;
  plot(he3_tab(p), p, 'b-');
  plot(he3_tc(p), p, 'b-');
  plot(he3_ta, he3_pa, 'bo');
  plot(he3_tb, he3_pb, 'bo');
  plot(he3_ts, he3_ps, 'bo');
  xlim([min(t) max(t)]);
  ylim([min(p) max(p)]);

  plot(he3_plts2gr(t), he3_pmelt_plts(t/1000), 'g-');

  xlabel('terature, mK');
  ylabel('pure, bar');

  print phase3.eps -deps "-S640,480" -color



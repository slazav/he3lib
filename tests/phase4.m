#!/usr/bin/octave -qf

  addpath ../matlab

  figure; clf; hold on;
  t=0.24:0.0001:0.51;

  plot(t, he3_pmelt_plts(t), 'r-');
%  plot(he3_tm_plts, he3_pm_plts, 'ro');

  plot(t, he3_pmelt(t), 'b-');
%  plot(he3_tm, he3_pm, 'bo');
  xlim([min(t) max(t)]);

  xlabel('terature, K');
  ylabel('pure, bar');

  print phase4.eps -deps "-S640,480" -color



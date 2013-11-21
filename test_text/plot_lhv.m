#!/usr/bin/octave -qf

  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;

  ttc = 0:0.01:1;
  p   = 0:5:30;

  for i=1:length(p)
    v = he3_text_lhv(ttc, p(i));
    plot(ttc, v, 'b.-');
  end

  plot(ttc, v, 'r.-');

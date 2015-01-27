#!/usr/bin/octave -qf

  addpath ../matlab
  figure; clf; hold on;

  ttc = 0:0.01:1;
  p   = 0:1:33;

  for i=1:length(p)
     v = he3_nu_b(ttc, p(i));
     plot(ttc, v, 'g.-');

     v1 = he3_nu_b1(ttc, p(i));
     plot(ttc, v1/2, 'b.-');
  end
  plot(ttc, v, 'r.-');



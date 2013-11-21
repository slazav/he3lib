#!/usr/bin/octave -qf

  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;
  ttc = 0:0.001:1;

  gap=he3_bcsgap(ttc);
  z3 = he3_z3(ttc, gap);
  z5 = he3_z5(ttc, gap);
  z7 = he3_z7(ttc, gap);

  ll = he3_lambda(ttc, gap);

  plot(ttc, z3, 'r-');
  plot(ttc, z5, 'g-');
  plot(ttc, z7, 'b-');
  plot(ttc, ll, 'm-');

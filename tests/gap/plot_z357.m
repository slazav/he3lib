#!/usr/bin/octave -qf

  addpath ../../matlab
  figure; clf; hold on;
  ttc = 0:0.001:1;

  gap=he3_bcsgap(ttc);
  z3 = he3_z3(ttc, gap);
  z5 = he3_z5(ttc, gap);
  z7 = he3_z7(ttc, gap);

  y0 = he3_yosida(ttc, gap, 0);
  y2 = he3_yosida(ttc, gap, 2);
  y5 = he3_yosida(ttc, gap, 5);
  y7 = he3_yosida(ttc, gap, 7);


  ll = he3_lambda(ttc, gap);

  plot(ttc, z3, 'r-');
  plot(ttc, z5, 'g-');
  plot(ttc, z7, 'b-');
%  plot(ttc, ll, 'm-');

  plot(ttc, (1-y0)./gap.^2, 'ro');
%  plot(ttc, y72, 'bo');

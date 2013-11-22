#!/usr/bin/octave -qf

  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;

  p   = 0:1:30;

  g = he3_gdk(p);     % K
  plot(p, g, 'b.-');

  g1 = he3_gd_exp(p); % 1e32 1/(erg cm^3)
  g1 = g1 .* he3_vm(p) / const_na / const_kb;
  plot(p, g1, 'r.-');

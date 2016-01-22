#!/usr/bin/octave -qf

  addpath ../../matlab
  figure; clf; hold on;
  ttc = 0:0.001:1;
  p=29;

  gap=he3_gap(ttc,p);
  y0 = he3_yosida(ttc, gap, 0);
  f0a = he3_f0a(p);

  chi1 = he3_chi_b(ttc,p);

  % full formula with f2a
  chi2 = @(f2a) (2/3 + (1/3+1/5*f2a)*y0)./...
         (1+f0a*(2/3+1/3*y0) + 1/5*f2a*(1/3+(2/3+f0a)*y0)).*...
         (1+f0a);

  plot(ttc, chi1, 'r-');
  plot(ttc, chi2(0), 'b-');
  plot(ttc, chi2(-0.051), 'r-');

  plot([0 1], [0.33 0.33], 'g-');



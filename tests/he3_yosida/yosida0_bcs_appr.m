#!/usr/bin/octave -qf

 % Einzel-2003 interpolation of yosida0 function (BCS)
 % Good for understanding low- and high-temperature behaviour.

 % for some reason accuracy worse then in the paper.
 % also formulas in the paper are not consistent with
 % the interpolation method.

  addpath ../../octave
  figure; clf; hold on;
  ttc = 0:0.005:1;
  gap=he3_bcsgap(ttc);
  y = he3_yosida(ttc, gap, 0);


  gap0 = pi/exp(const_euler) % BCS gap at T=0
  dcbcn0 = 12D0/7D0/const_z3

  % lowest-temperature approximation
  y00=sqrt(2*pi*gap0./ttc).*exp(-gap0./ttc);

  % next term in ttc/gap
  bet = 3.0/8.0;
  y0=y00.*(1 + bet*ttc./gap0);

  % value and derivative of y0 in ttc=1
  y00tc = sqrt(2*pi*gap0).*exp(-gap0);
  y0tc  = y00tc*(1 + bet./gap0);

  k= (s + 0.5 - gap0)/(1 - y0tc);
  yi = y0.*(1-ttc.^k) + exp(gap0-gap0./ttc) .* ttc.^(k-0.5);

  plot(ttc, y0./y, 'b-');
  plot(ttc, yi./y, 'r-');
  plot([0 1], [1 1]);

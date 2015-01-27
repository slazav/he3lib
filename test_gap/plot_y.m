#!/usr/bin/octave -qf

  addpath ../matlab
  figure; clf; hold on;
  ttc = 0:0.005:1;
  gap=he3_bcsgap(ttc);

  y = he3_yosida(ttc, gap, 0);

  plot(ttc, y, 'k-');

  y00=sqrt(2*pi*gap(1)./ttc).*exp(-gap(1)./ttc);
  y0=y00.*(1+3/8 * ttc./gap(1));

  k=2.388693;
  yi = y0.*(1-ttc.^k) + exp(gap(1)-gap(1)./ttc) .* ttc.^(k-1/2);

#  yi = y0.*(1-ttc.^k) + y0/y0(end) .* ttc.^k;

  plot(ttc, y00, 'r-');
  plot(ttc, y0, 'b-');
  plot(ttc, yi, 'g-');

  plot(ttc, yi./y, 'r-');
  plot([0 1], [1 1]);

#!/usr/bin/octave -qf
  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;

  ttc=0:0.01:1;
  p=0;
  h=200;

  dv=he3_text_ldv(ttc,p);
  hv=he3_text_lhv(ttc,p);

  plot(ttc, dv./(hv*h^2), 'b-')

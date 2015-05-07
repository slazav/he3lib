#!/usr/bin/octave -qf
  addpath ../matlab

  figure; clf; hold on;

  ttc=0.01;
  p=0:30;

  % recent theory
  c1=he3_cperp(ttc,p);
  c2=he3_cpar(ttc,p);

  % measured
  c3=(387./(20.544+p)+5.94)*100;
  c4=(348./(21.138+p)+5.13)*100;

  plot(p, c1, 'r-')
  plot(p, c2, 'b-')

  plot(p, c3, 'm-')
  plot(p, c4, 'c-')

#!/usr/bin/octave -qf
  addpath ../../matlab

  figure; clf; hold on;

  ttc=0.01:0.01:1;
  p=0;
  f1a = he3_f1a(p);
  chi = he3_chi_n(p) .* he3_chi_b(ttc,p);
  gap=he3_gap(ttc,p);
  y0=he3_yosida(ttc,gap,0);

 f3a=20;
  rl = he3_dorfle_rl(ttc,p, f1a, f3a);
  rt = he3_dorfle_rt(ttc,p, f1a, f3a);

  c =  2*he3_gyro^2./chi ...
        .* const_hbar^2./4/he3_amass./he3_meff(p) ...
        .* (1-y0).*he3_rho(p)/10;

  cper2 = c .* (2*rt+rl);
  cpar2 = c .* (4*rt);

  c1=he3_cperp(ttc,p);
  c2=he3_cpar(ttc,p);

  plot(ttc, c1, 'r-')
  plot(ttc, c2, 'b-')

  plot(ttc, sqrt(cper2), 'm-')
  plot(ttc, sqrt(cpar2), 'c-')

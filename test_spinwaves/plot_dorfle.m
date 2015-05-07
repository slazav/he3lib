#!/usr/bin/octave -qf
  addpath ../matlab

  figure; clf; hold on;

  ttc=0.01;
  p=0:30;
  f1a = he3_f1a(p);
  chi = he3_chi_n(p) .* he3_chi_b(ttc,p);
  y0=0;

  rl = he3_dorfle_rl(ttc,p, f1a, 0);
  rt = he3_dorfle_rt(ttc,p, f1a, 0);

  c =  2*he3_gyro^2./chi ...
        .* const_hbar^2./4/he3_amass./he3_meff(p) ...
        .* (1-y0).*he3_rho(p)/10;

  cper2 = c .* (2*rt+rl);
  cpar2 = c .* (4*rt);

  c1=he3_cperp(ttc,p);
  c2=he3_cpar(ttc,p);

  plot(p, c1, 'r-')
  plot(p, c2, 'b-')

  plot(p, sqrt(cper2), 'm-')
  plot(p, sqrt(cpar2), 'c-')

#!/usr/bin/octave -qf
  addpath ../../matlab

  figure; clf; hold on;

  ttc=0:0.01:1;
  p=0:5:30;

  cols = 'rgbcmkg';

  for i=1:length(p)
    a=he3_text_a(ttc,p(i));

    gape  = he3_gap(ttc,p(i)) * he3_tc(p(i)) * 1e-3 * const_kb;
    legge = const_hbar*he3_nu_b(ttc,p(i))*2*pi;


    chib = he3_chi_n(p(i)) * he3_chi_b(ttc, p(i));
    f0a  = he3_f0a(p(i));
    a1 = 5/36 * (legge./gape).^2 .*chib /(1+f0a)^2;

    plot(ttc, a, [cols(i) '-'])
    plot(ttc, a1, [cols(i) '.'])


  end

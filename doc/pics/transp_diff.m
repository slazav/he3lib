#!/usr/bin/octave -qf

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = [0.5:0.01:0.95 0.951:0.001:1];
  p=30;

  f=2e6;

  plot(ttc, he3_sdiff_hpar(ttc, p), 'r-');
  plot(ttc, he3_sdiff_hperp(ttc, p), 'r-', 'linewidth', 2);
  plot(ttc, he3_sdiff(ttc, p, f), 'b-', 'linewidth', 2);
  plot(ttc, he3_sdiff(ttc, p, 0), 'm--');

  ttc = 1:0.01:1.5;
  plot(ttc, he3_sdiff_nh(ttc, p), 'r-', 'linewidth', 2);
  plot(ttc, he3_sdiff_nperp(ttc, p, f), 'b-', 'linewidth', 2);

  plot([1 1], [0 0.05], 'k--');

  xlim([0.5 1.5])
  ylim([0 0.05])
  xlegend('T/T_c')
  ylegend('D, cm^2/s')

  text(0.95,0.010, 'D_\perp(1MHz)','fontweight','bold');
  text(1.20,0.021, 'D_\perp(0Hz) = \tau_{||}(0Hz)','fontweight','bold');
  text(1.20,0.019, 'D_{||}(1MHz)','fontweight','bold');
  text(0.66,0.040, 'D_{||}(0Hz)','fontweight','bold');
  text(0.66,0.038, 'D_{||}(1MHz)','fontweight','bold');
  text(0.90,0.040, 'D_\perp(0Hz)','fontweight','bold');

  % improved Samuli's code
  %addpath ~/he3lib/diff
  %s1d=diff_coeff(p, ttc, 1e4);
  %s2d=diff_coeff(p, ttc, 1e6);
%  s2=load('sdiff_sam2');
%  plot(s2.ttc, s2.s1d/dn, 'm-');
%  plot(s2.ttc, s2.s2d/dn, 'm-');

%   'Superfluid D_{perp} 10 kHz',...

  print -deps -color transp_diff.eps

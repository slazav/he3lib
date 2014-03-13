#!/usr/bin/octave -qf

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = [0.2:0.01:0.95 0.951:0.001:1 1.01:0.01 1.5];
  p=30;

  plot(ttc, he3_diff_hperp_zz(ttc, p), 'g-', 'linewidth', 1);
%  plot(ttc, he3_diff_hpar_zz(ttc, p), 'g-', 'linewidth', 2);

  for f=[0 1e6]
    plot(ttc, he3_diff_perp_zz(ttc, p, f), 'r-', 'linewidth', 2);
    plot(ttc, he3_diff_perp_xx(ttc, p, f), 'b-', 'linewidth', 2);
    plot(ttc, he3_diff_perp_zz_im(ttc, p, f), 'r-');
    plot(ttc, he3_diff_perp_xx_im(ttc, p, f), 'b-');
%    plot(ttc, he3_diff_par_zz(ttc, p, f), 'g-', 'linewidth', 2);
%    plot(ttc, he3_diff_par_xx(ttc, p, f), 'm-');
  end


  plot([1 1], [-1 1], 'k--');
  plot([0 2], [0 0], 'k--');

  xlim([0.2 1.5])
  ylim([-0.05 0.08])
  xlabel('T/T_c')
  ylabel('D, cm^2/s')

  text(1.17,0.025, 'D^\perp (0Hz)','fontweight','bold');
  text(1.10,0.008, 'D^\perp (1MHz)','fontweight','bold');
  text(0.45,0.067, 'D^\perp_{xx} (0Hz)','fontweight','bold');
  text(0.33,0.051, 'Re D^\perp_{xx} (1MHz)','fontweight','bold');
  text(0.74,0.066, 'D^\perp_{zz} (0Hz)','fontweight','bold');
  text(0.48,0.015, 'Re D^\perp_{zz} (1MHz)','fontweight','bold');
  text(0.46,-0.011, 'Im D^\perp_{xx} (0Hz)','fontweight','bold');
  text(0.73,-0.035, 'Im D^\perp_{zz} (1MHz)','fontweight','bold');

  % improved Samuli's code
  %addpath ~/he3lib/diff
  %s1d=diff_coeff(p, ttc, 1e4);
  %s2d=diff_coeff(p, ttc, 1e6);
%  s2=load('sdiff_sam2');
%  plot(s2.ttc, s2.s1d/dn, 'm-');
%  plot(s2.ttc, s2.s2d/dn, 'm-');

%   'Superfluid D_{perp} 10 kHz',...

  print transp_diff.eps -deps -color "-S800,600"

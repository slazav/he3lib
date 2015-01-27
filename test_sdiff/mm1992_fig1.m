function mm1992_fig1()
% Markelov, Mukharsky, PhB178 1992 fig1.
  addpath ../matlab
  figure; clf; hold on;
  ttc=0:0.01:1;

  plot(ttc, he3_diff_perp_xx(ttc, 30, 460000), 'r-', 'linewidth', 2);
  plot(ttc, he3_diff_perp_zz(ttc, 30, 460000), 'r-');
  plot(ttc, he3_diff_perp_xx_im(ttc, 30, 460000), 'b-', 'linewidth', 2);
  plot(ttc, he3_diff_perp_zz_im(ttc, 30, 460000), 'b-');
  plot(ttc, he3_diff_par_xx(ttc, 30, 460000), 'm-', 'linewidth', 2);
  plot(ttc, he3_diff_par_zz(ttc, 30, 460000), 'm-');
  plot([0 1], [0 0], 'k-');

  plotdat('mm1992_fig1_imperpxx.dat', 'k-', 'linewidth', 2);
  plotdat('mm1992_fig1_imperpzz.dat', 'k-');
  plotdat('mm1992_fig1_perpxx.dat', 'k-', 'linewidth', 2);
  plotdat('mm1992_fig1_perpzz.dat', 'k-');
  plotdat('mm1992_fig1_parxx.dat', 'k-', 'linewidth', 2);
  plotdat('mm1992_fig1_parzz.dat', 'k-');

  xlim([0 1.0]);
  ylim([-0.08 0.15]);
  xlabel('T/T_c');
  ylabel('D, cm^2/s');

  text(0.27, 0.030, 'Re D_{xx}^\perp');
  text(0.38, 0.015, 'Re D_{zz}^\perp');
  text(0.40,-0.013, 'Im D_{xx}^\perp');
  text(0.41,-0.047, 'Im D_{zz}^\perp');
  text(0.47, 0.100, 'D_{xx}^{||}');
  text(0.15, 0.081, 'D_{zz}^{||}');

  print -deps -color mm1992_fig1.eps

end

function plotdat(f, c)
  [x,y] = textread(f, '%f %f', 'commentstyle', 'shell');
  plot(x,y, c);
end


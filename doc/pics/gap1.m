#!/usr/bin/octave -qf

addpath ../../matlab

figure; clf; hold on;

  ttc=[0.001:0.001:1.1];

  gg00=he3_trivgap(ttc, 0);
  gg30=he3_trivgap(ttc, 30);

  plot(ttc, he3_bcsgap(ttc), 'k--')
  plot(ttc, gg00, 'k-')
  plot(ttc, gg30, 'k-', 'linewidth', 2)

  plot(ttc, he3_yosida(ttc,gg00,0), 'r')
  plot(ttc, he3_yosida(ttc,gg00,2), 'g')
  plot(ttc, he3_yosida(ttc,gg00,4), 'b')

  plot(ttc, he3_yosida(ttc,gg30,0), 'r', 'linewidth', 2)
  plot(ttc, he3_yosida(ttc,gg30,2), 'g', 'linewidth', 2)
  plot(ttc, he3_yosida(ttc,gg30,4), 'b', 'linewidth', 2)

  plot(ttc, he3_chi_b(ttc,0), 'k');
  plot(ttc, he3_chi_b(ttc,30), 'k', 'linewidth', 2);

  plot(ttc, he3_rho_nb(ttc,0), 'b');
  plot(ttc, he3_rho_nb(ttc,30), 'b', 'linewidth', 2);

  legend('BCS gap', 'P=0', 'P=30', 'location', 'northeast')

  text(0.8, 1.5, '\Delta', 'fontsize', 8, 'color', 'k');
  text(0.51, 0.28, 'Y_0', 'fontsize', 8, 'color', 'r');
  text(0.68, 0.20, 'Y_2', 'fontsize', 8, 'color', 'g');
  text(0.87, 0.12, 'Y_4', 'fontsize', 8, 'color', 'b');
  text(0.55, 0.84, '\rho^n_B/\rho_N', 'fontsize', 8, 'color', 'b');
  text(0.23, 0.45, '\chi_B/\chi_N'   , 'fontsize', 8, 'color', 'k');
  xlim([0 1]);
  ylim([0 2]);
  set(gca,'yTick', 0:0.2:2);
  grid on;
  print gap1.eps -deps "-S640,480" -color


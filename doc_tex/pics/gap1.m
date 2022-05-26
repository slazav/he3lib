#!/usr/bin/octave-cli -qf

  pkg load he3lib
  figure; clf; hold on;
  subplot(1,2,1); hold on;
  xlabel "T/Tc"
  title "Gap/Tc"

  ttc=[0.001:0.001:1.1];

  gg00=he3_trivgap(ttc, 0);
  gg30=he3_trivgap(ttc, 30);

  plot(ttc, he3_bcsgap(ttc), 'k--')
  plot(ttc, gg00, 'k-')
  plot(ttc, gg30, 'k-', 'linewidth', 2)

  xlim([0 1]);
  ylim([0 2]);
  set(gca, 'yTick', 0:0.2:2);
  legend('bcsgap', 'trivgap P=0', 'trivgap P=30', 'location', 'southwest')

  grid on;

  subplot(1,2,2); hold on;
  xlabel "T/Tc"
  title "Gap-related functions"

  plot(ttc, he3_yosida(ttc,gg00,0), 'r')
  plot(ttc, he3_yosida(ttc,gg00,2), 'color', [0 0.6 0] )
  plot(ttc, he3_yosida(ttc,gg00,4), 'b')

  plot(ttc, he3_yosida(ttc,gg30,0), 'r', 'linewidth', 2)
  plot(ttc, he3_yosida(ttc,gg30,2), 'color', [0 0.6 0], 'linewidth', 2)
  plot(ttc, he3_yosida(ttc,gg30,4), 'b', 'linewidth', 2)

  plot(ttc, he3_z3(ttc,gg00), 'm', 'linewidth', 1)
  plot(ttc, he3_z5(ttc,gg00), 'color', [0.6 0.6 0], 'linewidth', 1)
  plot(ttc, he3_z7(ttc,gg00), 'color', [0 0.6 0.6], 'linewidth', 1)

  plot(ttc, he3_z3(ttc,gg30), 'm', 'linewidth', 2)
  plot(ttc, he3_z5(ttc,gg30), 'color', [0.6 0.6 0], 'linewidth', 2)
  plot(ttc, he3_z7(ttc,gg30), 'color', [0 0.6 0.6], 'linewidth', 2)


  plot(ttc, he3_chi_b(ttc,0), 'k');
  plot(ttc, he3_chi_b(ttc,30), 'k', 'linewidth', 2);

  plot(ttc, he3_rho_nb(ttc,0), 'b');
  plot(ttc, he3_rho_nb(ttc,30), 'b', 'linewidth', 2);


  text(0.8, 1.5, '\Delta', 'fontsize', 8, 'color', 'k');
  text(0.70, 0.35, 'Y_0', 'fontsize', 8, 'color', 'r');
  text(0.75, 0.25, 'Y_2', 'fontsize', 8, 'color', [0 0.6 0]);
  text(0.87, 0.12, 'Y_4', 'fontsize', 8, 'color', 'b');
  text(0.70, 0.92, '\rho^n_B/\rho', 'fontsize', 8, 'color', 'b');
  text(0.10, 0.42, '\chi_B/\chi_N'   , 'fontsize', 8, 'color', 'k');

  text(0.10, 0.95, 'Z_3'   , 'fontsize', 8, 'color', 'm');
  text(0.10, 0.72, 'Z_5'   , 'fontsize', 8, 'color', [0.6 0.6 0]);
  text(0.10, 0.58, 'Z_7'   , 'fontsize', 8, 'color', [0 0.6 0.6]);

  xlim([0 1]);
  ylim([0 1]);
  set(gca, 'yTick', 0:0.1:1);
  grid on;
  print gap1.eps -deps "-S500,180" -color


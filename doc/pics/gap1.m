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

  legend('BCS gap', 'P=0', 'P=30', 'location', 'northeast')

%  plot(ttc, he3_z3(ttc,gg), 'r')
%  plot(ttc, he3_z5(ttc,gg), 'g')
%  plot(ttc, he3_z7(ttc,gg), 'b')
  text(0.8, 1.5, '\Delta', 'fontsize', 8, 'color', 'k');
  text(0.73, 0.7,  'Y_0', 'fontsize', 8, 'color', 'r');
  text(0.83, 0.5, 'Y_2', 'fontsize', 8, 'color', 'g');
  text(0.93, 0.2, 'Y_4', 'fontsize', 8, 'color', 'b');
  xlim([0 1.02]);


print gap1.eps -deps "-S640,480" -color


function plot_yosida_lt()
  addpath ../matlab
  figure; clf; hold on;
  ttc = 0.01:0.01:1;
  p=1;

  gap=he3_trivgap(ttc,p);
  Z3 = he3_z3(ttc,gap);
  Z5 = he3_z5(ttc,gap);
  Z7 = he3_z7(ttc,gap);

  plot(ttc, 1-he3_yosida(ttc,gap, 0), 'g-', 'linewidth', 2);
  plot(ttc, Z3, 'g-');
  plot(ttc, Z5, 'b-');
  plot(ttc, Z7, 'r-');

  xlim([0 1]);
  ylim([0 2]);
  xlabel('ttc')
  print -deps -color plot_z.eps
end

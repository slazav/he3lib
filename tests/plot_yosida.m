function plot_yosida_lt()
  t='Yosida functions Y_n';
  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;
  ttc = 0.01:0.01:1;
  p=1;

  gap=he3_trivgap(ttc,p);

  Y0 = he3_yosida(ttc,gap,0);
  Y1 = he3_yosida(ttc,gap,1);
  Y2 = he3_yosida(ttc,gap,2);

  plot(ttc, Y0, 'r-');
  plot(ttc, Y1, 'g-');
  plot(ttc, Y2, 'b-');

  xlabel('ttc')
  ylabel('Y_n')
  title(t);
  legend('Y0', 'Y1', 'Y3');
  print -deps -color plot_yosida.eps
end

function einzel1991_fig5()
% Einzel JLTP84 (1991) p.343 fig.5
  t='The generalized Yosida functions Y_n(T)';
  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;
  ttc = 0.01:0.01:1;
  p=1;

  gap=he3_trivgap(ttc,p);

  Y0 = he3_yosida(ttc,gap,0);
  Y1 = he3_yosida(ttc,gap,1);
  Y2 = he3_yosida(ttc,gap,2);
  Y3 = he3_yosida(ttc,gap,3);
  Y4 = he3_yosida(ttc,gap,4);

  plot(ttc, Y0, 'r-');
  plot(ttc, Y1, 'g-');
  plot(ttc, Y2, 'b-');
  plot(ttc, Y3, 'm-');
  plot(ttc, Y4, 'c-');

  legend('n=0', '1', '2', '3', '4',...
    'location', 'northwest');
  xlim([0.2 1.0]);
  ylim([0 1.0]);
  xlabel('T/T_c');
  ylabel('Y_n(T)');
  title(t);
  print -deps -color einzel1991_fig5.eps
end

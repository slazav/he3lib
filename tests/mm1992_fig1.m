function mm1992_fig1()
  % Einzel JLTP84 (1991) p.353 fig.7

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.1:0.005:1;
  f=460000;

  plot(ttc, he3_sdiff_mm(ttc, 30, f), 'r--');
  plot(ttc, he3_sdiff(ttc, 30, f), 'b--');

  legend('Markelov-Mukharsky', 'Bunkov-Einzel')

%  ylim([0 0.5]);
  xlabel('T/T_c');
  ylabel('D^\perp, cm^2/s');

  print -deps -color mm1992_fig1.eps
end
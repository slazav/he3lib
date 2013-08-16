function einzel_plot5()
  % Einzel JLTP84 (1991) p.353 fig.7

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.4:0.01:1;
  f=460000;

  plot(ttc, he3_sdiff(ttc, 29, f), 'r-');
  plot(ttc, he3_sdiff(ttc, 20, f), 'g-');
  plot(ttc, he3_sdiff(ttc, 11, f), 'b-');
  plot(ttc, he3_sdiff(ttc,  0, f), 'm-');

  legend(...
    '26 bar',...
    '20 bar',...
    '11 bar',...
    ' 0 bar',...
    ''
  );
  ylim([0 0.5]);
  xlabel('T/T_c');
  ylabel('D, cm^2/s');

  print -deps -color einzel_plot5.eps
end
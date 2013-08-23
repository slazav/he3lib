function einzel1991_fig8()
  % Einzel JLTP84 (1991) p.353 fig.7

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.4:0.01:1;
  f=460000;

  plot(ttc, he3_sdiff(ttc, 29, f), 'r-');
  plot(ttc, he3_sdiff(ttc, 20, f), 'g-');
  plot(ttc, he3_sdiff(ttc, 11, f), 'b-');
  plot(ttc, he3_sdiff(ttc,  0, f), 'm-');

  [x,y] = textread('einzel1991_fig8c1.dat', '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'k-');
  [x,y] = textread('einzel1991_fig8c2.dat', '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'k-');
  [x,y] = textread('einzel1991_fig8c3.dat', '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'k-');
  [x,y] = textread('einzel1991_fig8c4.dat', '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'k-');

  [x,y] = textread('einzel1991_fig8e1.dat', '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'ko');
  [x,y] = textread('einzel1991_fig8e2.dat', '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'kd');
  [x,y] = textread('einzel1991_fig8e3.dat', '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'ko');
  [x,y] = textread('einzel1991_fig8e4.dat', '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'ks');

  legend(...
    '29 bar',...
    '20 bar',...
    '11 bar',...
    ' 0 bar',...
    ''
  );
  xlim([0.4 1]);
  ylim([0 0.5]);
  xlabel('T/T_c');
  ylabel('D, cm^2/s');

  print -deps -color einzel1991_fig8.eps
end
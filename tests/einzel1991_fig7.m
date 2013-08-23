function einzel1991_fig7()
  % Einzel JLTP84 (1991) p.353 fig.7

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0:0.01:1;
  p=30;
  dn = he3_sdiff(1, p, 0);

  [x,y] = textread('einzel1991_fig7a.dat',...
    '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'k-');
  [x,y] = textread('einzel1991_fig7b.dat',...
    '%f %f', 'commentstyle', 'shell');
  plot(x,y, 'k-');

  plot(ttc, he3_sdiff(ttc, p, 1e4)/dn, 'r-');
  plot(ttc, he3_sdiff(ttc, p, 1e6)/dn, 'b-');
  plot(ttc, he3_sdiff_hperp(ttc, p)/dn, 'g-');
  plot(ttc, he3_sdiff(ttc, p, 0)/dn, 'm-');

  legend(...
    'original plot, 10kHz',...
    'original plot, 1000kHz',...
    '10 kHz',...
    '1000 kHz',...
    'Hydrodynamic',...
    '0Hz',...
    ''
  );
  ylim([0 10]);
  xlabel('T/T_c');
  ylabel('D/D0');

  print -deps -color einzel1991_fig7.eps
end
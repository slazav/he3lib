function einzel1991_fig7()
  % Einzel JLTP84 (1991) p.353 fig.7

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0:0.01:1;
  p=30;

  dn = he3_sdiff(1, p, 0);
  plot(ttc, he3_sdiff_hperp(ttc, p)/dn, 'k-');
  plot(ttc, he3_sdiff(ttc, p, 0)/dn, 'k-');
  plot(ttc, he3_sdiff(ttc, p, 1e4)/dn, 'r-');
  plot(ttc, he3_sdiff(ttc, p, 1e6)/dn, 'b-');

  legend(...
    'Hydrodynamic',...
    '0Hz',...
    '10 kHz',...
    '1000 kHz',...
    ''
  );
  ylim([0 10]);
  xlabel('T/T_c');
  ylabel('D/D0');

  print -deps -color einzel1991_fig7.eps
end
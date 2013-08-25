function einzel1991_fig8()
  % Einzel JLTP84 (1991) p.353 fig.7

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.4:0.01:1;
  f=460000;

  plot(ttc, he3_diff(ttc, 29, f), 'r-');
  plot(ttc, he3_diff(ttc, 20, f), 'g-');
  plot(ttc, he3_diff(ttc, 11, f), 'b-');
  plot(ttc, he3_diff(ttc,  0, f), 'm-');

  plotdat('einzel1991_fig8c1.dat', 'k-');
  plotdat('einzel1991_fig8c2.dat', 'k-');
  plotdat('einzel1991_fig8c3.dat', 'k-');
  plotdat('einzel1991_fig8c4.dat', 'k-');
  plotdat('einzel1991_fig8e1.dat', 'ko');
  plotdat('einzel1991_fig8e2.dat', 'kd');
  plotdat('einzel1991_fig8e3.dat', 'ko');
  plotdat('einzel1991_fig8e4.dat', 'ks');

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

function plotdat(f, c)
  [x,y] = textread(f, '%f %f', 'commentstyle', 'shell');
  plot(x,y, c);
end

#!/usr/bin/octave-cli -qf

  % collect and fit all experimental data
  pkg load he3lib
  figure; clf; hold on;

  subplot(1,2,1); hold on;
    P=[0 10 20 30];
    f0=500:10:900;
    Imin=2;
    col='rgbm';
    leg={};
    for i=1:length(P)
      plot(f0, rota_qball_fz0(P(i),f0*1000,Imin), [ col(i) '-' ], 'linewidth', 2)
      leg{end+1} = sprintf('P=%d bar', P(i));
    end
    for i=1:length(P)
      plot(f0, rota_qball_fr0(P(i),f0*1000,Imin), [ col(i) '--' ], 'linewidth', 2)
    end
    legend(leg, 'location', 'northeast')
    xlabel('f_0, kHz');
    title('f_r^0, f_z^0, Hz');
    text(700,350,'f_r^0')
    text(550,180,'f_z^0')

  subplot(1,2,2); hold on;
    P=[0 10 20 30];
    f0=500:10:900;
    Imin=2;
    col='rgbm';
    leg={};
    for i=1:length(P)
      plot(f0, rota_qball_az0(P(i),f0*1000,Imin), [ col(i) '-' ], 'linewidth', 2)
      leg{i} = sprintf('P=%d bar', P(i));
    end
    for i=1:length(P)
      plot(f0, rota_qball_ar0(P(i),f0*1000,Imin), [ col(i) '--' ], 'linewidth', 2)
      leg{i} = sprintf('P=%d bar', P(i));
    end
    legend(leg, 'location', 'northeast')
    xlabel('f0, kHz');
    title('a_r^0,a_z^0, cm');
    text(600,0.055,'a_z^0')
    text(550,0.030,'a_r^0')

  print qball.eps -deps "-S500,250" -color



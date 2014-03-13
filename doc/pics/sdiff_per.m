#!/usr/bin/octave -qf

  addpath ../../matlab
  figure; clf; hold on;

  p=0;
  f=[1 2 5];
  freqs=[f*0.1 f f*10 f*1e2 f*1e3 f*1e4 f*1e5 f*1e6];
  ttc=0.1:0.01:1.2;
  semilogy(ttc, he3_diff_hperp_zz(ttc, p), 'g-', 'linewidth', 2);
  for f=freqs;
    semilogy(ttc, he3_diff_perp_zz(ttc, p, f), 'b-');
    semilogy(ttc, he3_diff_perp_xx(ttc, p, f), 'r-');
    if f>5e4 && f<=1e6
      text(1.05, he3_diff_perp_xx(1.05, p, f), sprintf('f=%.1e',f))
    end
  end

  legend('he3\_diff\_hperp\_zz', 'he3\_diff\_perp\_zz', 'he3\_diff\_perp\_xx');
  title('D_{perp} vs T/T_c at P=0')

  xlim([0.1 1.2]);
  ylim([0.03 30]);

  xlabel('T/T_c');
  ylabel('D_\perp');

  print -deps -color sdiff_per.eps

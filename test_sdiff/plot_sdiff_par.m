function plot_sdiff_par()
  % Einzel JLTP84 (1991) p.353 fig.7

  addpath ../matlab

  figure; clf; hold on;

  p=0;
  f=[1 2 5];
  freqs=[f*0.1 f f*10 f*1e2 f*1e3 f*1e4 f*1e5 f*1e6];
  ttc=0.1:0.01:1.2;
  semilogy(ttc, he3_diff_hperp_zz(ttc, p), 'g-', 'linewidth', 2);
  for f=freqs;
    semilogy(ttc, he3_diff_par_zz(ttc, p, f), 'r-');
    semilogy(ttc, he3_diff_par_xx(ttc, p, f), 'b-');
  end
  legend('he3\_diff\_hpar\_zz', 'he3_diff\_par\_zz', 'he3_diff\_par\_xx');

  xlim([0.1 1.2]);
  ylim([1 30]);

  xlabel('T/T_c');
  ylabel('D_\par');

  print -deps -color plot_sdiff_par.eps
end

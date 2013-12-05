function plot_sdiff_par()
  % Einzel JLTP84 (1991) p.353 fig.7

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;

  p=0;
  o=[1 2 5];
  om=[o*0.1 o o*10 o*1e2 o*1e3 o*1e4 o*1e5 o*1e6];
  ttc=0.1:0.01:1.2;
  for o=om;
    semilogy(ttc, he3_diff_par_zz(ttc, p, o), 'r-');
    semilogy(ttc, he3_diff_par_xx(ttc, p, o), 'b-');
  end
  semilogy(ttc, he3_diff_hpar_zz(ttc, p), 'g-');

  xlim([0.1 1.2]);
  ylim([1 30]);

  xlabel('T/T_c');
  ylabel('D_\perp');

  print -deps -color plot_sdiff_par.eps
end
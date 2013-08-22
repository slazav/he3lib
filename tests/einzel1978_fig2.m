function einzel1978_fig2()
  % Einzel JLTP32 (1978) p.36 fig.2
  t = 'Mean free path of Bogoliubov quasiparticle in He3-B';

  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;
  ttc = 0:0.005:1;
  p=0;

  tn0=he3_tau_n0(1, p);
  vf=he3_vf(p);
  gap=he3_trivgap(ttc,p);
  l=he3_fpath(ttc, p);

  plot(1-ttc, l.*exp(-gap./ttc)/vf/tn0 , 'k.-');

  % Low temperature approximation
  g0 = he3_scatt_g0(p);
  d0 = he3_scatt_d0(p);
  w0 = (1 - 2/3*g0 + d0);
  ttch = 0:0.01:0.2;
  gaph=he3_trivgap(ttch,p);
  l_lt = sqrt(2*pi)/3./gaph.^2/w0;
  plot(1-ttch, l_lt, 'r-');

%  plot(ttc, l, 'k-');
%  plot(ttc, l_lt, 'b');

  xlim([0 1]);
  ylim([0 0.8]);
  xlabel('1-T/T_c')
  ylabel('L exp(-\Delta/T)/(v_F \tau_N(0,T_c))')
  title(t);

  print -deps -color einzel1978_fig2.eps

end

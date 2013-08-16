#!/usr/bin/octave -qf

  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;
  ttc = 0:0.005:1;
  p=0;

  tn0=he3_tau_n0(1, p);
  vf=he3_vf(p);
  gap=he3_trivgap(ttc,p);

  l=he3_fpath(ttc, p);

  sc=exp(-gap./ttc)/vf/tn0;

  % Low temperature approximation
  w0 = he3_scatt_w0(p);
  l_lt = sqrt(2*pi)/3./gap.^2/ w0;


  plot(1-ttc, sc.*l, 'k.-');
  plot(1-ttc, l_lt, 'r-');

%  plot(ttc, l, 'k-');
%  plot(ttc, l_lt, 'b');

  xlim([0 1]);
  ylim([0 0.8]);
%  xlabel('T/T_c')
%  ylabel('L exp(-\Delta/T)/(v_F \tau_N(0,T_c))')
  title('Mean free path of Bogoliubov quasiparticle in He3-B')

  print -deps -color plot_fpath.eps

function einzel1991_fig6()
  % Einzel JLTP84 (1991) p.351 fig.6
  % Hydrodynamic spin diffusion
  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0:0.001:1;
  p=0;

  tn0=he3_tau_n0(1, p);
  vf=he3_vf(p);
  gap=he3_trivgap(ttc,p);
  l=he3_fpath(ttc, p);


%  dtc = he3_vf(p)^2/3/he3_chi_b(1,p) * he3_tau_n_av(1,p)
  dtc = he3_sdiff_hpar(1,p);
  Dpar  = he3_sdiff_hpar(ttc,p);
  Dperp = he3_sdiff_hperp(ttc,p);

  plot(ttc, Dpar/dtc, 'b.-');
  plot(ttc, Dperp/dtc, 'r.-');

  xlim([0 1])
  ylim([0 24])

  print -deps -color einzel1991_fig6a.eps

  figure; clf; hold on;
  plot(ttc, Dpar/dtc, 'b.-');
  plot(ttc, Dperp/dtc, 'r.-');
  xlim([0.96 1])
  ylim([0.8 1.2])

  print -deps -color einzel1991_fig6b.eps

end

function einzel1991_fig6()
  % Einzel JLTP84 (1991) p.351 fig.6
  % Hydrodynamic spin diffusion
  addpath ../../matlab

  figure; clf; hold on;
  ttc = 0:0.001:1;
  p=0;

  tn0=he3_tau_n0(1, p);
  vf=he3_vf(p);
  gap=he3_trivgap(ttc,p);
  l=he3_fpath(ttc, p);


  dtc = he3_diff_hpar_zz(1,p);
  Dpar  = he3_diff_hpar_zz(ttc,p);
  Dperp = he3_diff_hperp_zz(ttc,p);

  plotdat('einzel1991_fig6a1.dat', 'k-');
  plotdat('einzel1991_fig6a2.dat', 'k-');
  plot(ttc, Dpar/dtc, 'b.-');
  plot(ttc, Dperp/dtc, 'r.-');

  xlim([0 1])
  ylim([0 24])

  print -deps -color einzel1991_fig6a.eps

  figure; clf; hold on;
  plotdat('einzel1991_fig6b1.dat', 'k-');
  plotdat('einzel1991_fig6b2.dat', 'k-');
  plot(ttc, Dpar/dtc, 'b.-');
  plot(ttc, Dperp/dtc, 'r.-');
  xlim([0.96 1])
  ylim([0.8 1.2])

  print -deps -color einzel1991_fig6b.eps

end
function plotdat(f, c)
  [x,y] = textread(f, '%f %f', 'commentstyle', 'shell');
  plot(x,y, c);
end

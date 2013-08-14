function einzel_plot2()
  % Einzel JLTP32 (1978) p.36 fig.2
  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.00:0.01:1;
  p=0;

  tn0=he3_tau_n0(1, p);
  vf=he3_vf(p);
  gap=he3_trivgap(ttc,p);
  l=he3_fpath(ttc, p);
%  l = vf*tn0*(he3_yosida(ttc,gap,2)./he3_yosida(ttc,gap,0)).^2;

  g0 = he3_scatt_g0(p);
  d0 = he3_scatt_d0(p);
  w0 = (1 - 2/3*g0 + d0);
  l_lt = vf*tn0*sqrt(2*pi)/3./gap.^2 .*exp(gap./ttc)/w0;

  plot(1-ttc, l.*exp(-gap./ttc)/vf/tn0, 'k-');
%  plot(1-ttc, l_lt.*exp(-gap./ttc)/vf/tn0, 'r-');

  plot(1-ttc, sqrt(2*pi)/3./gap.^2/w0, 'r-');

  plot(1-ttc, he3_yosida0(ttc, gap).*sqrt(ttc).*exp(gap./ttc), 'g-');

%  plot(ttc, l, 'k-');
%  plot(ttc, l_lt, 'b');

  xlim([0 1]);
  ylim([0 10]);

end

function einzel_plot1()
  % Einzel JLTP32 (1978) p.35 fig.1
  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.00:0.01:1;
  p=0;

  tn=he3_tau_n0(ttc,p);
  gap=he3_trivgap(ttc,p);
  t0=he3_tau0(ttc, p);
  ta=he3_tau_av(ttc, p);

  plot(gap./ttc, tn./t0, 'k-');
  plot(gap./ttc, tn./ta/6, 'r-');
  xlabel('\Delta/T')
  ylabel('\tau_N(0,T)/\tau')
  legend('Lifetime at the Fermi level',...
         'Thermal average lifetime * 6')

  g0 = he3_scatt_g0(p);
  d0 = he3_scatt_d0(p);
  w0 = (1 - 2/3*g0 + d0);
  cq = (3*d0-g0)/3.0/w0 - 3/4.0;

  t0lt = 1/3/w0 * sqrt(2*pi*(gap./ttc).^3) .* exp(gap./ttc) .* tn .*...
    (1 + cq*ttc./gap);

  plot(gap./ttc, tn./t0lt, 'b-')

  ylim([0 1]);
  xlim([0 9]);


  print -deps -color einzel_plot1.eps

end

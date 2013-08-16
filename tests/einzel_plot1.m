function einzel_plot1()
  % Einzel JLTP32 (1978) p.35 fig.1
  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.00:0.01:1;
  p=0;

  gap=he3_trivgap(ttc,p);

  tn0=he3_tau_n0(ttc,p);
  tna=he3_tau_n_av(ttc,p);
  t0=he3_tau0(ttc, p);
  ta=he3_tau_av(ttc, p);

  plot(gap./ttc, tn0./t0, 'b-');
  plot(gap./ttc, tna./ta, 'r-');
  xlabel('\Delta/T')
  ylabel('\tau_N/\tau')
  legend('\tau_N(0)/\tau(0)',...
         '<\tau_N(E)>/<\tau(E)>'...
  )

  ylim([0 1]);
  xlim([0 9]);

  print -deps -color einzel_plot1.eps

end

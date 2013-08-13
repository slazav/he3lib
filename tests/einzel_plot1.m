function plot_tauav()
  % Einzel JLTP32 (1978) p.35 fig.1
  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.00:0.01:1;
  p=0;

  tn=he3_tau_n0tc(p) ./ ttc.^2;
  gap=he3_trivgap(ttc,p);
  t0=he3_tau0(ttc, p);

  plot(gap./ttc, tn./t0, 'k-');
  xlim([0 9]);

end
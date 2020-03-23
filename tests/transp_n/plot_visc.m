function plot_visc()

  find_figure('spin diff'); clf; hold on;

  p=0;
  ttc=0.8:0.01:1.3;

  T = 1e-3*ttc*he3_tc(p); % T, K

  % Weatley-1976 \eta T^2 [poise*mK^2] vs P [bar]
  pp  = [0 3 6 9 12 15 18 21 24 27 30 33 34.36];
  et2 = [1.834 1.73 1.63 1.54 1.46 1.38 1.30 1.22 1.14 1.06 0.99 0.92 0.88];
  plot(pp, et2, "r*")

  p = 0:0.2:33;
  ttc = 1;
  et2a = he3_hvisc(ttc, p).*(ttc*he3_tc(p)).^2;
  plot(p, et2a, "b-");
end




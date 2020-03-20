function plot_sd()

  find_figure('spin diff'); clf; hold on;

  p=30;
  ttc=0:0.01:1;

  t = 1e-3*ttc*he3_tc(p); % T, K

  % experiment, whwatley 1978,  see VW p.36
  TE = 0.3e-6 ./ (t*1000).^2;

  % in he3lib
  T0 = 4/3 * he3_tau_n0(ttc,p);

  a0s = he3_f0s(p)/(1+he3_f0s(p))
  a1s = he3_f1s(p)/(1+he3_f1s(p)/3.0)
  a0a = he3_f0a(p)/(1+he3_f0a(p))
  a1a = he3_f1a(p)/(1+he3_f1a(p)/3.0)

  W = abs(a0a)^2 + 3*abs(a1a)^2

  semilogy(ttc, TE, "b-")
  semilogy(ttc, T0, "r-")

end




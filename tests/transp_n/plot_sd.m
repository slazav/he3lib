function plot_sd()

  find_figure('spin diff'); clf; hold on;

  p=0;
  ttc=0.8:0.01:1.3;

  T = 1e-3*ttc*he3_tc(p); % T, K

  % gas-kynetic viscosity, heat capacity, spin-diffusion coefficients
  % VW 2.40

  tn0 = pi/4 * he3_tau_n0(ttc,p);

  E0 = 1/5.0 * he3_rho(p) * he3_vf(p)^2 * tn0;
  K0 = 1/3.0 * he3_gammaf(p) * he3_vf(p)^2 * tn0 .* T;
  D0 = 1/3.0 * he3_vf(p)^2 * tn0;

  D1 = he3_diffn_hydr(ttc,p);

  la1p = he3_scatt_l1a(p);
  ls1m = 1;
  ls2p = 0;

  % VW 2.69, 2.71
  E1a = E0 * he3_mm(p) * f_eta(ls2p);
  K1a = K0 * f_k(ls1m);
  D1a = D0 * (1 + he3_f0a(p)) * f_eta(la1p);

  l=0:0.1:2;
  semilogy(ttc, E0, "b--")
  semilogy(ttc, K0, "r--")
  semilogy(ttc, D0, "g--")

  semilogy(ttc, E1a, "b-")
  semilogy(ttc, K1a, "r-")
  semilogy(ttc, D1a, "g-")

  semilogy(ttc, D1, "g.")
%  semilogy(ttc, D1a, "g-")

end




function e1990_cmp()

  % complare table in Einzel-1990 paper with he3lib

  pp = [0 10 20 30];
  f1s = [5.39 8.57 11.09 13.5];
  f2s = [-1.26 0.215 0.461 0.712];
  vf  = [5976 4630 3920 3429];
  tau0 = [0.500 0.120 0.054 0.040]  * 1e-6;
  dccn = [1.43 1.70 1.83 1.91];
  gap0 = [1.76 1.82 2.00 2.10];
  l2  = [0.68 0.74 0.7 0.74];
  g0  = [0.1 0.1 0.1 0.1];
  d0  = [0.3 0.3 0.3 0.3];

  find_figure('e1990_cmp'); clf; hold on;

  plot(pp, f1s./he3_f1s(pp), 'r*-')
  plot(pp, f2s./he3_f2s(pp), 'b*-')  % ref.29
  plot(pp, vf./he3_vf(pp), 'bo-')
  plot(pp, dccn./he3_dcbcn(pp), 'b.-')
  plot(pp, gap0./he3_gap(0, pp), 'g.-')
%  plot(pp, tau0./he3_tau0(1, pp), 'mo-')
%  plot(pp, g0./he3_scatt_g0(pp), 'mx-')
%  plot(pp, d0./he3_scatt_d0(pp), 'rx-')
end
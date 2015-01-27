function plot_yosida_lt()
  t='Check low temperature behaviour of Y0';
  addpath ../matlab
  figure; clf; hold on;
  ttc = 0.01:0.005:1;
  p=1;

  gap=he3_trivgap(ttc,p);
  Y0  = he3_yosida(ttc,gap,0);
  Y0o = 1-he3_z3(ttc,gap);
  Y0lt = sqrt(2*pi*gap./ttc) .* exp(-gap./ttc) .* (1+0.375*ttc./gap(1));

  Y2  = he3_yosida(ttc,gap,2);
  Y2o = 1 - he3_z5(ttc,gap);
  Y2lt = ttc./gap .* sqrt(2*pi*gap./ttc).*exp(-gap./ttc) .* (1-ttc./gap(1));

  Y4  = he3_yosida(ttc,gap,4);
  Y4o = 1 - he3_z7(ttc,gap);
  Y4lt = 3*(ttc./gap).^2 .* sqrt(2*pi*gap./ttc).*exp(-gap./ttc) .* (1-16/3*ttc./gap(1));

%  addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Diffusion coeff/Calc_diff_coeff');
%  for i=1:length(ttc); Y0s(i) = yosida(p, ttc(i), 0); end

  plot([0 0.1], [1 1], 'k-');
%  plot(ttc, Y0o./Y0lt, 'g-');
%  plot(ttc, Y0./Y0lt, 'r-');
%  plot(ttc, Y0s./Ylt, 'b-');

%  plot(ttc, Y2o./Y2lt, 'g-');
%  plot(ttc, Y2./Y2lt, 'r-');

%  plot(ttc, Y4o./Y4lt, 'g-');
%  plot(ttc, Y4./Y4lt, 'r-');

  plot(ttc, 1-he3_z3(ttc,gap), 'g-');
  plot(ttc, he3_z5(ttc,gap), 'b-');
  plot(ttc, he3_z7(ttc,gap), 'c-');
  plot(ttc, Y2, 'r-');

  xlim([0 1]);
  ylim([0 1]);
  xlabel('ttc')
  ylabel('Y0/Y0lt')
  title(t);
  legend('LT limit', 'he3lib', 'texture lib', 'Samuli''s code')
  print -deps -color plot_yosida_lt.eps
end

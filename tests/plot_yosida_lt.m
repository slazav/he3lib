function plot_yosida_lt()
  t='Check low temperature behaviour of Y0';
  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;
  ttc = 0.01:0.01:1;
  p=1;

  gap=he3_trivgap(ttc,p);
  Y0  = he3_yosida(ttc,gap,0);
  Y0o = he3_yosida0(ttc,gap);
  Ylt = sqrt(2*pi*gap./ttc) .* exp(-gap./ttc);

%  addpath('/rota/Analysis/NMRcalc/spinwaves/Spinwave_relaxation/Diffusion coeff/Calc_diff_coeff');
%  for i=1:length(ttc); Y0s(i) = yosida(p, ttc(i),0); end

  plot([0 0.1], [1 1], 'k-');
  plot(ttc, Y0o./Ylt, 'g-');
  plot(ttc, Y0./Ylt, 'r-');
  plot(ttc, he3_yosida0_fast(ttc, gap)./Ylt, 'b-');

%  plot(ttc, Y0, 'r-');
%  plot(ttc, he3_yosida0_fast(ttc, gap), 'b-');
%  plot(ttc, Y0s./Ylt, 'b-');

  xlim([0 1]);
  ylim([0 1.2]);
  xlabel('ttc')
  ylabel('Y0/Y0lt')
  title(t);
  legend('correct value', 'he3lib', 'texture lib', 'Samuli''s code')
  print -deps -color plot_yosida_lt.eps
end

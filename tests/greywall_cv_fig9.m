function greywall_cv_fig9()
% Greywall PRB2 1983 fig9.
  addpath ~/he3lib/lib/matlab
  figure; clf; hold on;
  t=0.16:0.01:2.5;

  plot(t, he3_cv_n(t, 36.82), 'r-');
  plot(t, he3_cv_n(t, 32.59), 'r-');
  plot(t, he3_cv_n(t, 30.39), 'r-');
  plot(t, he3_cv_n(t, 28.89), 'r-');
  plot(t, he3_cv_n(t, 27.70), 'r-');
  plot(t, he3_cv_n(t, 26.84), 'r-');
  plot(t, he3_cv_n(t, 26.17), 'r-');
  xlim([0 2.5]);
  ylim([0.2 0.9]);
  xlabel('T, K');
  ylabel('C_v/R');

  print -deps -color greywall_cv_fig9.eps

end

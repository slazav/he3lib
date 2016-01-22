function greywall_cv_fig14()
% Greywall PRB2 1983 fig14.
  addpath ../../matlab
  figure; clf; hold on;
  t=(0:0.01:50)/1000;

  plot(t*1000, he3_cv_n(t, he3_vm(32.50))./t, 'r-');
  plot(t*1000, he3_cv_n(t, he3_vm(29.30))./t, 'r-');
  plot(t*1000, he3_cv_n(t, he3_vm(22.22))./t, 'r-');
  plot(t*1000, he3_cv_n(t, he3_vm(17.01))./t, 'r-');
  plot(t*1000, he3_cv_n(t, he3_vm(11.00))./t, 'r-');
  plot(t*1000, he3_cv_n(t, he3_vm(5.040))./t, 'r-');
  plot(t*1000, he3_cv_n(t, he3_vm(0.060))./t, 'r-');

  xlim([0 50]);
  ylim([2.5 4.5]);
  xlabel('T, mK');
  ylabel('C_v/RT');
  set(gca, 'yTick', [2.4:0.2:4.4]);

  print -deps -color greywall_cv_fig14.eps

end

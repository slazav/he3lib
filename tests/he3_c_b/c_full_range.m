function full_range()

  addpath ../../octave
  find_figure('full_range'); clf; hold on;

  for P=[0 10 20 30]
    tc = he3_tc(P)
    t=(tc:0.01:10)/1000;  % 0..10mK
    ttc = 0.05:0.01:1;

    tt = [ttc*tc t*1000];
    cc = [he3_c_b(ttc, P) he3_cv_n(t, he3_vm(P))];

    plot(tt, cc, 'r-');
    plot(t*1000, he3_gammaf(P)*t, 'r--');
  end

  xlabel('T, mK');
  ylabel('C_v/R');
  set(gca, 'yTick', [2.4:0.2:4.4]);

  %print -deps -color greywall_cv_fig14.eps

end

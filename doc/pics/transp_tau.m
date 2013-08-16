#!/usr/bin/octave -qf

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.00:0.01:1;
  p=0;

  tn=he3_tau_n0(ttc, p);
  Y0=he3_yosida(ttc, p, 0);

  if 0
    % low temp
    tav_lt=he3_tau_av(ttc, p);
    save tau_lt tav_lt
  else
    load tau_lt
  end

  if 0
    % high temp
    tav_ht=he3_tau_av(ttc, p);
    save tau_ht tav_ht
  else
    load tau_ht
  end

  ttcs = 0.00:0.1:1;
  if 0
    % original Samuli's code
    addpath('s')
    for i = 1:length(ttc)
      ts(i) = tau(p, ttc(i));
    end
    save tau_samuli ts
  else
    load tau_samuli
  end

  semilogy(ttc, Y0.*he3_tau_n0(ttc, p), 'r-');
  semilogy(ttc, Y0.*he3_tau_n_av(ttc, p), 'b-');

  semilogy(ttc, Y0.*he3_tau0(ttc, p), 'r-', 'linewidth', 2);
  semilogy(ttc, Y0.*he3_tau_av(ttc, p), 'b-', 'linewidth', 2);

  semilogy(ttc, Y0.*he3_tau0lt(ttc, p), 'r-.');

  semilogy(ttc, Y0.*tav_lt, 'b--');
  semilogy(ttc, Y0.*tav_ht, 'b-.');

  semilogy(ttc, Y0.*ts, 'g-', 'linewidth', 2);

  xlim([0.3 1]);
  ylim(10.^[-6.5 -5]);

  legend(...
    '\tau_N(0)',...
    '<\tau_N(E)>',...
    '\tau(0)',...
    '<\tau(E)>',...
    '\tau(E) low temp limit',...
    '<\tau(E)> with low temp I(E,T)',...
    '<\tau(E)> with high temp I(E,T)',...
    'Samuli''s code',...
    ''
  );
  xlabel('T/T_c');
  ylabel('Y_0 \tau');
  print -deps -color transp_tau.eps

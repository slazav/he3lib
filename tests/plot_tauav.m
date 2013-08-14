function plot_tauav()
  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0.00:0.01:1;
  p=0;

  tn=he3_tau_n0(ttc, p);

  if 0
    % full range temp
    tav=he3_tau_av(ttc, p);
    save tau_av tav
  else
    load tau_av
  end

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
    save ts ts
  else
    load ts
  end

  semilogy(ttc, tav./tn, 'r-');
  semilogy(ttc, tav_lt./tn, 'b-');
  semilogy(ttc, tav_ht./tn, 'g-');
  semilogy(ttc, ts./tn, 'm--');
  semilogy(ttc, he3_tau0(ttc, p)./tn, 'k-');
  semilogy([0 1], [1 1], 'k-');
  ylim(10.^[-1 2]);

  legend(...
    '\tau',...
    'low temp approx',...
    'high temp approx',...
    'Samuli''s code',...
    '\tau_0',...
    '\tau_N'...
  );
  xlabel('T/T_c');
  ylabel('\tau/\tau_N');
  print -deps -color plot_tauav.eps
end
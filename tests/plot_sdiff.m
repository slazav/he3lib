function plot_sdiff()
  addpath ~/he3lib/lib/matlab
  addpath ~/he3lib/diff

  figure; clf; hold on;
  ttc = 0:0.01:1;
  p=30;

  dn = he3_sdiff(1, p, 0);
  plot(ttc, he3_sdiff(ttc, p, 1e4)/dn, 'r-');
  plot(ttc, he3_sdiff(ttc, p, 1e6)/dn, 'b-');

  s1d=diff_coeff(p, ttc, 1e4);
  s2d=diff_coeff(p, ttc, 1e6);
  plot(ttc, s1d/dn, 'r-');
  plot(ttc, s2d/dn, 'b-');
  save sdiff s1d s2d ttc p


%  plot(ttc, he3_sdiff(ttc, 29, 460000), 'b-');

%  ylim(10.^[-1 2]);
%
%  legend(...
%    '\tau',...
%    'low temp approx',...
%    'high temp approx',...
%    'Samuli''s code',...
%    '\tau_0',...
%    '\tau_N'...
%  );
%  xlabel('T/T_c');
%  ylabel('\tau/\tau_N');
%  print -deps -color plot_sdiff.eps
end
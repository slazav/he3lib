function plot_sdiff()
  addpath ~/he3lib/lib/matlab
  addpath ~/he3lib/diff

  figure; clf; hold on;
  ttc = 0:0.01:1;
  p=0;
  no0=823000;

  plot(ttc, he3_sdiff(ttc, p, 1e4), 'r-');
  plot(ttc, he3_sdiff(ttc, p, 1e6), 'b-');
%  plot(ttc, diff_coeff(p, ttc, no0), 'b-');

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
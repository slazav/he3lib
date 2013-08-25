#!/usr/bin/octave -qf

  addpath ~/he3lib/lib/matlab

  figure; clf; hold on;
  ttc = 0:0.01:1;

  semilogy(ttc, he3_fpath(ttc, 0),  'r-');
  semilogy(ttc, he3_fpath(ttc, 10), 'g-');
  semilogy(ttc, he3_fpath(ttc, 20), 'b-');
  semilogy(ttc, he3_fpath(ttc, 30), 'm-');

  xlim([0 1])
  ylim(10.^[-5 2])
  xlabel('T/T_c')
  ylabel('mean free path, cm')
  legend('0 bar', '10 bar', '20 bar', '30 bar');
%  set(gca, 'xTick', [0:0.1:1]);
  grid on;

  print transp_fp.eps -deps -color "-S640,480"


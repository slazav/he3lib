#!/usr/bin/octave -qf

  % collect and fit all experimental data

  addpath ../../matlab

function do_fit(p, ttc, nub, c)
  xx=0.001:0.001:1;
  plot(xx, he3_nu_b(xx, p).^2, [c(1) '-']);
  plot(ttc,  nub, c)
end


  figure; clf; hold on;
  xlabel('T/T_c');
  ylabel('f_B^2, Hz^2');

  % data from Hakonen, et al., LJTP 76, p.225-283 (1989) 
  [ttc005 nub005] = textread('nu_b_data/1989-Hakonen-005.dat', '%f %f', 'commentstyle', 'shell');
  [ttc050 nub050] = textread('nu_b_data/1989-Hakonen-050.dat', '%f %f', 'commentstyle', 'shell');
  [ttc102 nub102] = textread('nu_b_data/1989-Hakonen-102.dat', '%f %f', 'commentstyle', 'shell');
  [ttc155 nub155] = textread('nu_b_data/1989-Hakonen-155.dat', '%f %f', 'commentstyle', 'shell');
  [ttc250 nub250] = textread('nu_b_data/1989-Hakonen-250.dat', '%f %f', 'commentstyle', 'shell');
  do_fit(0.5, ttc005, nub005, 'r*');
  do_fit(5.0, ttc050, nub050, 'g*');
  do_fit(10.2, ttc102, nub102, 'b*');
  do_fit(15.5, ttc155, nub155, 'm*');
  do_fit(25.0, ttc250, nub250, 'c*');

  % Ahonen 1976 (Wheatley-78)
  [ttc187 nub187] = textread('nu_b_data/1976-Ahonen-187.dat', '%f %f', 'commentstyle', 'shell');
  [ttc211 nub211] = textread('nu_b_data/1976-Ahonen-211.dat', '%f %f', 'commentstyle', 'shell');
  [ttc254 nub254] = textread('nu_b_data/1976-Ahonen-254.dat', '%f %f', 'commentstyle', 'shell');
  [ttc290 nub290] = textread('nu_b_data/1976-Ahonen-290.dat', '%f %f', 'commentstyle', 'shell');
  [ttc320 nub320] = textread('nu_b_data/1976-Ahonen-320.dat', '%f %f', 'commentstyle', 'shell');
  do_fit(18.7, ttc187, nub187, 'ro');
  do_fit(21.1, ttc211, nub211, 'go');
  do_fit(25.4, ttc254, nub254, 'bo');
  do_fit(29.0, ttc290, nub290, 'mo');
  do_fit(32.0, ttc320, nub320, 'co');

  %  ROTA measurements (vs MCT) 5-6/08/2003
  [ttc102n nub102n] = textread('nu_b_data/2003-Rota-102.dat',  '%f %f', 'commentstyle', 'shell');
  nub102n = nub102n*1e10;
  %  ROTA measurements 2007, counterflow peak
  [ttc290a nub290a] = textread('nu_b_data/2007-Rota-290.dat',  '%f %f', 'commentstyle', 'shell');
  nub290a = nub290a.^2;
  do_fit(10.2, ttc102n, nub102n, 'g.');
  do_fit(29.0, ttc290a, nub290a, 'rv');

  %Rota 2011-2013
  [ttc005l nub005l] = textread('nu_b_data/2011-Rota-005.dat', '%f %f', 'commentstyle', 'shell');
  [ttc082 nub082]   = textread('nu_b_data/2011-Rota-082.dat', '%f %f', 'commentstyle', 'shell');
  [ttc156 nub156]   = textread('nu_b_data/2011-Rota-156.dat', '%f %f', 'commentstyle', 'shell');
  ttc156 = ttc156(13:end); nub156 = nub156(13:end);
  do_fit(0.5, ttc005l, nub005l, 'rd');
  do_fit(8.2, ttc082, nub082, 'bd');
  do_fit(15.6, ttc156, nub156, 'gd');

  print nu_b.eps -deps "-S800,600" -color



#!/usr/bin/octave-cli -qf

  % collect and fit all experimental data

  pkg load he3lib

function [press, y] = do_fit(press, ttc, nub, c)
  nubt = ft1(ttc, press);
  y = sum(nub.*nubt)/sum(nubt.*nubt); % fit nub = y * ft1(ttc) via y
  xx=0.001:0.001:1;
  subplot(3,2,1:4);
  plot(xx, y*ft1(xx, press), [c(1) '-']);
%  plot(xx, ft1(xx, press), [c(1) '--']);
  plot(ttc,  nub, c)
  text(0.02 + (c(1)-95)/200, y*ft1(0, press)+1e9, num2str(press));
  subplot(3,2,5:6);
  plot(press, y, c);
end

% theoretical function
function y = ft1(ttc, p)
  y = he3_nu_b1(ttc, p).^2;
end


  figure; clf;
  subplot(3,2,1:4); hold on;
  xlabel('T/T_c');
  ylabel('f_B^2, Hz^2');
  subplot(3,2,5:6); hold on;
  xlabel('pressure, bar');

  % data from Hakonen, et al., LJTP 76, p.225-283 (1989) 
  [ttc005 nub005] = textread('nu_b_data/1989-Hakonen-005.dat', '%f %f', 'commentstyle', 'shell');
  [ttc050 nub050] = textread('nu_b_data/1989-Hakonen-050.dat', '%f %f', 'commentstyle', 'shell');
  [ttc102 nub102] = textread('nu_b_data/1989-Hakonen-102.dat', '%f %f', 'commentstyle', 'shell');
  [ttc155 nub155] = textread('nu_b_data/1989-Hakonen-155.dat', '%f %f', 'commentstyle', 'shell');
  [ttc250 nub250] = textread('nu_b_data/1989-Hakonen-250.dat', '%f %f', 'commentstyle', 'shell');
  [P(end+1) A(end+1)] = do_fit(0.5, ttc005, nub005, 'r*');
  [P(end+1) A(end+1)] = do_fit(5.0, ttc050, nub050, 'g*');
  [P(end+1) A(end+1)] = do_fit(10.2, ttc102, nub102, 'b*');
  [P(end+1) A(end+1)] = do_fit(15.5, ttc155, nub155, 'm*');
  [P(end+1) A(end+1)] = do_fit(25.0, ttc250, nub250, 'c*');

  % Ahonen 1976 (Wheatley-78)
  [ttc187 nub187] = textread('nu_b_data/1976-Ahonen-187.dat', '%f %f', 'commentstyle', 'shell');
  [ttc211 nub211] = textread('nu_b_data/1976-Ahonen-211.dat', '%f %f', 'commentstyle', 'shell');
  [ttc254 nub254] = textread('nu_b_data/1976-Ahonen-254.dat', '%f %f', 'commentstyle', 'shell');
  [ttc290 nub290] = textread('nu_b_data/1976-Ahonen-290.dat', '%f %f', 'commentstyle', 'shell');
  [ttc320 nub320] = textread('nu_b_data/1976-Ahonen-320.dat', '%f %f', 'commentstyle', 'shell');
  [P(end+1) A(end+1)] = do_fit(18.7, ttc187, nub187, 'ro');
  [P(end+1) A(end+1)] = do_fit(21.1, ttc211, nub211, 'go');
  [P(end+1) A(end+1)] = do_fit(25.4, ttc254, nub254, 'bo');
  [P(end+1) A(end+1)] = do_fit(29.0, ttc290, nub290, 'mo');
  [P(end+1) A(end+1)] = do_fit(32.0, ttc320, nub320, 'co');

  %  ROTA measurements (vs MCT) 5-6/08/2003
  [ttc102n nub102n] = textread('nu_b_data/2003-Rota-102.dat',  '%f %f', 'commentstyle', 'shell');
  nub102n = nub102n*1e10;
  %  ROTA measurements 2007, counterflow peak
  [ttc290a nub290a] = textread('nu_b_data/2007-Rota-290.dat',  '%f %f', 'commentstyle', 'shell');
  nub290a = nub290a.^2;
  [P(end+1) A(end+1)] = do_fit(10.2, ttc102n, nub102n, 'g.');
  [P(end+1) A(end+1)] = do_fit(29.0, ttc290a, nub290a, 'rv');

  %Rota 2011-2013
  [ttc005l nub005l] = textread('nu_b_data/2011-Rota-005.dat', '%f %f', 'commentstyle', 'shell');
  [ttc082 nub082]   = textread('nu_b_data/2011-Rota-082.dat', '%f %f', 'commentstyle', 'shell');
  [ttc156 nub156]   = textread('nu_b_data/2011-Rota-156.dat', '%f %f', 'commentstyle', 'shell');
  ttc156 = ttc156(13:end); nub156 = nub156(13:end);
  [P(end+1) A(end+1)] = do_fit(0.5, ttc005l, nub005l, 'rd');
  [P(end+1) A(end+1)] = do_fit(8.2, ttc082, nub082, 'bd');
  [P(end+1) A(end+1)] = do_fit(15.6, ttc156, nub156, 'gd');

  subplot(3,2,1:4); hold on;
%  xlim([0.8 1]);
%  ylim([0 3e10]);

  print nu_b.eps -deps "-S1200,1600" -color



#!/usr/bin/octave-cli -qf

  pkg load he3lib

  # copper parameters:
  Bint = 0.36e-3;    # [T], dipolar feld in copper
  gyro = 71.118e6;   # [rad/s/T] gyromagnetic ratio of copper
  spin = 1.5;        # spin 3/2

  t=exp([log(0.001):0.01:log(20)]);     # mK

  B0 = 0;    # T
  B1 = 0.08; # T
  B2 = 0.8;  # T
  B3 = 8;    # T

  x0 = gyro*const_hbar*sqrt(B0^2+Bint^2) ./ (const_kb * 1e-3 * t);
  x1 = gyro*const_hbar*sqrt(B1^2+Bint^2) ./ (const_kb * 1e-3 * t);
  x2 = gyro*const_hbar*sqrt(B2^2+Bint^2) ./ (const_kb * 1e-3 * t);
  x3 = gyro*const_hbar*sqrt(B3^2+Bint^2) ./ (const_kb * 1e-3 * t);

  figure; clf;
  subplot(1,3,1); hold on;
  title("Magnetization")
  xlabel('T, mK');
  ylabel('M/Ms');

  semilogx(t,magn_par_m(x0, spin), 'k-');
  semilogx(t,magn_par_m(x1, spin), 'm-');
  semilogx(t,magn_par_m(x2, spin), 'b-');
  semilogx(t,magn_par_m(x3, spin), 'r-');
  xlim([0,t(end)])
  legend('B=0', 'B=0.08T', 'B=0.8T', 'B=8T', 'Location', 'SouthWest')


  subplot(1,3,2); hold on;
  title("Entropy")
  xlabel('T, mK');
  ylabel('S/R');

  semilogx(t,magn_par_s(x0, spin), 'k-');
  semilogx(t,magn_par_s(x1, spin), 'm-');
  semilogx(t,magn_par_s(x2, spin), 'b-');
  semilogx(t,magn_par_s(x3, spin), 'r-');
  xlim([0,t(end)])
  ylim([0,1.5])


  subplot(1,3,3); hold on;
  title("Heat capacity")
  xlabel('T,mK');
  ylabel('C/R');

  semilogx(t,magn_par_c(x0, spin), 'k-');
  semilogx(t,magn_par_c(x1, spin), 'm-');
  semilogx(t,magn_par_c(x2, spin), 'b-');
  semilogx(t,magn_par_c(x3, spin), 'r-');
  xlim([0,t(end)])
  ylim([0,0.8])

  print magn_par.png -dpng "-S800,300"



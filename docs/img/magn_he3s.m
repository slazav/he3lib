#!/usr/bin/octave-cli -qf

  pkg load he3lib

  # copper parameters:
  Bint = 700e-3     # [T], dipolar feld in solid helium-3
  gyro = 203.789e6  # [rad/s/T] gyromagnetic ratio, helium-3
  spin = 0.5        # spin 1/2, helium-3

  tc = 0.5 # mK

  t=exp([log(0.1):0.01:log(50)]);     # mK

  B0 = 0;     # T
  B1 = 0.08;  # T
  B2 = 0.8;   # T
  B3 = 8;     # T

  # arguments for paramagnetic functions:
  x0 = gyro*const_hbar*sqrt(B0^2+Bint^2) ./ (const_kb * 1e-3 * t);
  x1 = gyro*const_hbar*sqrt(B1^2+Bint^2) ./ (const_kb * 1e-3 * t);
  x2 = gyro*const_hbar*sqrt(B2^2+Bint^2) ./ (const_kb * 1e-3 * t);
  x3 = gyro*const_hbar*sqrt(B3^2+Bint^2) ./ (const_kb * 1e-3 * t);

  # arguments for CW functions
  y0 = spin * gyro*const_hbar*B0 ./ (const_kb * 1e-3 * tc);
  y1 = spin * gyro*const_hbar*B1 ./ (const_kb * 1e-3 * tc);
  y2 = spin * gyro*const_hbar*B2 ./ (const_kb * 1e-3 * tc);
  y3 = spin * gyro*const_hbar*B3 ./ (const_kb * 1e-3 * tc);


  figure; clf;
  subplot(1,3,1); hold on;
  title("Magnetization")
  xlabel('T, mK');
  ylabel('M/Ms');

  semilogx(t,magn_cw_m(t/tc,y0), 'k--');
  semilogx(t,magn_cw_m(t/tc,y1), 'm--');
  semilogx(t,magn_cw_m(t/tc,y2), 'b--');
  semilogx(t,magn_cw_m(t/tc,y3), 'r--');

  semilogx(t,magn_par_m(x0, spin), 'k-');
  semilogx(t,magn_par_m(x1, spin), 'm-');
  semilogx(t,magn_par_m(x2, spin), 'b-');
  semilogx(t,magn_par_m(x3, spin), 'r-');
  xlim([0,t(end)])
  ylim([0,1])


  subplot(1,3,2); hold on;
  title("Entropy")
  xlabel('T, mK');
  ylabel('S/R');

  semilogx(t,magn_cw_s(t/tc,y0), 'k--');
  semilogx(t,magn_cw_s(t/tc,y1), 'm--');
  semilogx(t,magn_cw_s(t/tc,y2), 'b--');
  semilogx(t,magn_cw_s(t/tc,y3), 'r--');

  semilogx(t,magn_par_s(x0, spin), 'k-');
  semilogx(t,magn_par_s(x1, spin), 'm-');
  semilogx(t,magn_par_s(x2, spin), 'b-');
  semilogx(t,magn_par_s(x3, spin), 'r-');
  xlim([0,t(end)])
  ylim([0,0.75])


  subplot(1,3,3); hold on;
  title("Heat capacity")
  xlabel('T,mK');
  ylabel('C/R');

  semilogx(t,magn_cw_c(t/tc,y0), 'k--');
  semilogx(t,magn_cw_c(t/tc,y1), 'm--');
  semilogx(t,magn_cw_c(t/tc,y2), 'b--');
  semilogx(t,magn_cw_c(t/tc,y3), 'r--');

  semilogx(t,magn_par_c(x0, spin), 'k-');
  semilogx(t,magn_par_c(x1, spin), 'm-');
  semilogx(t,magn_par_c(x2, spin), 'b-');
  semilogx(t,magn_par_c(x3, spin), 'r-');
  xlim([0,t(end)])
  ylim([0,1.5])
  legend('B=0T', 'B=0.08T', 'B=0.8T', 'B=8T', 'Location', 'NorthEast')

  print magn_he3s.png -dpng "-S800,300"




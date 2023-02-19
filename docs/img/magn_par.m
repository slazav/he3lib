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

  figure; clf;
  subplot(1,3,1); hold on;
  title("Magnetization")
  xlabel('T, mK');
  ylabel('M/Ms');

  Ms = magn_par_m(0, B3, Bint, gyro, spin);

  semilogx(t,magn_par_m(1e-3*t, B0, Bint, gyro, spin)/Ms, 'k-');
  semilogx(t,magn_par_m(1e-3*t, B1, Bint, gyro, spin)/Ms, 'm-');
  semilogx(t,magn_par_m(1e-3*t, B2, Bint, gyro, spin)/Ms, 'b-');
  semilogx(t,magn_par_m(1e-3*t, B3, Bint, gyro, spin)/Ms, 'r-');
  xlim([0,t(end)])
  legend('B=0', 'B=0.08T', 'B=0.8T', 'B=8T',...
         'Location', 'West')


  subplot(1,3,2); hold on;
  title("Entropy")
  xlabel('T, mK');
  ylabel('S/R');

  semilogx(t,magn_par_s(1e-3*t, B0, Bint, gyro, spin), 'k-');
  semilogx(t,magn_par_s(1e-3*t, B1, Bint, gyro, spin), 'm-');
  semilogx(t,magn_par_s(1e-3*t, B2, Bint, gyro, spin), 'b-');
  semilogx(t,magn_par_s(1e-3*t, B3, Bint, gyro, spin), 'r-');
  xlim([0,t(end)])
  ylim([0,1.5])


  subplot(1,3,3); hold on;
  title("Heat capacity")
  xlabel('T,mK');
  ylabel('C/R');

  semilogx(t,magn_par_c(1e-3*t, B0, Bint, gyro, spin), 'k-');
  semilogx(t,magn_par_c(1e-3*t, B1, Bint, gyro, spin), 'm-');
  semilogx(t,magn_par_c(1e-3*t, B2, Bint, gyro, spin), 'b-');
  semilogx(t,magn_par_c(1e-3*t, B3, Bint, gyro, spin), 'r-');
  xlim([0,t(end)])
  ylim([0,0.8])

  print magn_par.png -dpng "-S800,300" "-F:6"


available_graphics_toolkits
#!/usr/bin/octave-cli -qf

  pkg load he3lib


  t=exp(log(0.1):0.01:log(10))*1e-3;

  figure; clf;
  subplot(2,2,1); hold on;
  title("Magnetization")
  xlabel('T, mK');
  ylabel('M');

  Tc = 0.5e-3;
  B0 = 0;
  B1 = 0.08;
  B2 = 0.8;
  B3 = 8;
  gyro = 203.789e6

  semilogx(t*1e3,magn_cw_m(t,B0,Tc,gyro), 'k-');
  semilogx(t*1e3,magn_cw_m(t,B1,Tc,gyro), 'm-');
  semilogx(t*1e3,magn_cw_m(t,B2,Tc,gyro), 'b-');
  semilogx(t*1e3,magn_cw_m(t,B3,Tc,gyro), 'r-');
  ylim([0,70000])



  subplot(2,2,2); hold on;
  title("Susceptibility")
  xlabel('T,mK');
  ylabel('Chi');

  semilogx(t*1e3,magn_cw_chi(t,B0,Tc,gyro), 'k-');
  semilogx(t*1e3,magn_cw_chi(t,B1,Tc,gyro), 'm-');
  semilogx(t*1e3,magn_cw_chi(t,B2,Tc,gyro), 'b-');
  semilogx(t*1e3,magn_cw_chi(t,B3,Tc,gyro), 'r-');
  ylim([0,200000])
  legend('B=0', 'B=80mT', 'B=0.8T', 'B=8T', 'Location', 'NorthEast')



  subplot(2,2,3); hold on;
  title("Entropy")
  xlabel('T,mK');
  ylabel('S/R');

  semilogx(t*1e3,magn_cw_s(t,B0,Tc,gyro), 'k-');
  semilogx(t*1e3,magn_cw_s(t,B1,Tc,gyro), 'm-');
  semilogx(t*1e3,magn_cw_s(t,B2,Tc,gyro), 'b-');
  semilogx(t*1e3,magn_cw_s(t,B3,Tc,gyro), 'r-');
  ylim([0,0.8])



  subplot(2,2,4); hold on;
  title("Heat capacity")
  xlabel('T,mK');
  ylabel('C/R');

  semilogx(t*1e3,magn_cw_c(t,B0,Tc,gyro), 'k-');
  semilogx(t*1e3,magn_cw_c(t,B1,Tc,gyro), 'm-');
  semilogx(t*1e3,magn_cw_c(t,B2,Tc,gyro), 'b-');
  semilogx(t*1e3,magn_cw_c(t,B3,Tc,gyro), 'r-');
  ylim([0,1.5])

  print magn_cw.png -dpng "-S800,600" "-F:6"


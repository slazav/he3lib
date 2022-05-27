#!/usr/bin/octave-cli -qf

  pkg load he3lib


  t=0:0.005:1.7;

  figure; clf;
  subplot(2,2,1); hold on;
  title("Magnetization")
  xlabel('T/Tc');
  ylabel('M/mu');

  plot(t,magn_cw_m(t,0), 'k-');
  plot(t,magn_cw_m(t,0.01), 'm-');
  plot(t,magn_cw_m(t,0.1), 'b-');
  plot(t,magn_cw_m(t,1), 'r-');
  xlim([0,t(end)])
  legend('muB/kTc=0', 'muB/kTc=0.01', 'muB/kTc=0.1', 'muB/kTc=1', 'Location', 'SouthWest')


  subplot(2,2,2); hold on;
  title("Susceptibility")
  xlabel('T/Tc');
  ylabel('Chi Tc/mu^2');

  plot(t,magn_cw_chi(t,0), 'k-');
  plot(t,magn_cw_chi(t,0.01), 'm-');
  plot(t,magn_cw_chi(t,0.1), 'b-');
  plot(t,magn_cw_chi(t,1), 'r-');
  ylim([0,20])
  xlim([0,t(end)])



  subplot(2,2,3); hold on;
  title("Entropy")
  xlabel('T/Tc');
  ylabel('S/R');

  plot(t,magn_cw_s(t,0), 'k-');
  plot(t,magn_cw_s(t,0.01), 'm-');
  plot(t,magn_cw_s(t,0.1), 'b-');
  plot(t,magn_cw_s(t,1), 'r-');
  xlim([0,t(end)])



  subplot(2,2,4); hold on;
  title("Heat capacity")
  xlabel('T/Tc');
  ylabel('C/R');

  plot(t,magn_cw_c(t,0), 'k-');
  plot(t,magn_cw_c(t,0.01), 'm-');
  plot(t,magn_cw_c(t,0.1), 'b-');
  plot(t,magn_cw_c(t,1), 'r-');
  xlim([0,t(end)])


  print magn_cw.png -dpng "-S800,600"



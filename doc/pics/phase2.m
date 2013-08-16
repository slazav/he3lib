#!/usr/bin/octave -qf

  addpath ../../matlab

  figure; clf; hold on;
  temp=10.^(-3:0.1:1.5);
  press=10.^(-6:0.1:4);
  loglog(temp, he3_pmelt(temp), 'b-');
  loglog(he3_tc(press)/1e3, press, 'b-');
  loglog(he3_tab(press)/1e3, press, 'b-');
  loglog(he3_ta/1e3, he3_pa, 'r.');
  loglog(he3_tb/1e3, he3_pb, 'r.');
  loglog(he3_tabn/1e3, he3_pabn, 'r.');
  loglog(he3_tneel/1e3, he3_pneel, 'g.');
  loglog(temp, he3_pvap(temp), 'b-');
  loglog(he3_tsmin, he3_psmin, 'ro');
  loglog(he3_tcr, he3_pcr, 'ro');

  text(he3_tsmin, 10, 'he3\_psmin,');
  text(he3_tsmin,  3, 'he3\_tsmin');
  text(he3_tcr*1.2, 0.10, 'he3\_pcr,');
  text(he3_tcr*1.2, 0.03, 'he3\_tcr');
  text(0.005,80, 'he3\_pmelt(t)');
  text(0.4,1e-4,  'he3\_pvap(t)');

  text(0.01,0.01, 'Liquid','fontweight','bold');
  text(0.01,1000, 'Solid','fontweight','bold');
  text(5,0.003, 'Gas','fontweight','bold');
  text(0.0011,10, 'B','fontweight','bold');

  xlim([min(temp) max(temp)]);
  ylim([10^-5 10^4]);
  xlabel('temperature, K');
  ylabel('pressure, bar');

  print phase2.eps -deps "-S640,480" -color


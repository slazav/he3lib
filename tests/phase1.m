#!/usr/bin/octave -qf

  addpath ../matlab

  figure; clf; hold on;
  temp=0.3:0.01:3.5;
  plot(temp, he3_pmelt(temp/1000), 'b-');
  press=0:0.1:36;
  plot(he3_tab(press), press, 'b-');
  plot(he3_tc(press), press, 'b-');
  plot(he3_tabn, he3_pabn, 'ro');
  plot(he3_ta, he3_pa, 'ro');
  plot(he3_tb, he3_pb, 'ro');
  plot(he3_ts, he3_ps, 'go');
  xlim([0.5 3]);
  ylim([0 45]);
  text(he3_ts-0.2, he3_ps+5, 'he3\_ps,');
  text(he3_ts-0.2, he3_ps+2.5, 'he3\_ts');
  text(he3_tb-0.2, he3_pb+5, 'he3\_pb,');
  text(he3_tb-0.2, he3_pb+2.5, 'he3\_tb');
  text(he3_ta, he3_pa+5, 'he3\_pa,');
  text(he3_ta, he3_pa+2.5, 'he3\_ta');
  text(he3_tabn+0.1, he3_pabn, 'he3\_pabn,');
  text(he3_tabn+0.1, he3_pabn-2.5, 'he3\_tabn');
  text(1.2,10,   'he3\_tc(p)');
  text(1.5,26, 'he3\_tab(p)');

  text(2.17,30, 'A','fontweight','bold');
  text(1.3,20, 'B','fontweight','bold');
  text(2,5, 'Normal','fontweight','bold');
  text(2.1,42, 'Solid','fontweight','bold');
  xlabel('temperature, mK');
  ylabel('pressure, bar');

  print phase1.eps -deps "-S640,480" -color



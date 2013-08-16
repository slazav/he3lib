#!/usr/bin/octave -qf


function title1(t)
  text(0.5,0.85, t, ...
    'units', 'normalized', 'horizontalalignment', 'center');
end

  addpath ../../matlab

  figure; clf; hold on;
  temp=0.3:0.01:3.5;
  plot(temp, he3_pmelt(temp/1000), 'b-');
  press=0:0.1:36;
  plot(he3_tab(press), press, 'b-');
  plot(he3_tc(press), press, 'b-');
  plot(he3_tabn, he3_pabn, 'ro');
  plot(he3_ta, he3_pa, 'ro');
  plot(he3_tb, he3_pb, 'ro');
  plot(he3_tneel, he3_pneel, 'go');
  xlim([0.5 3]);
  ylim([0 45]);
  text(he3_tneel-0.2, he3_pneel+3, 'he3\_pneel,');
  text(he3_tneel-0.2, he3_pneel+1.5, 'he3\_tneel');
  text(he3_tb-0.2, he3_pb+3, 'he3\_pb,');
  text(he3_tb-0.2, he3_pb+1.5, 'he3\_tb');
  text(he3_ta, he3_pa+3, 'he3\_pa,');
  text(he3_ta, he3_pa+1.5, 'he3\_ta');
  text(he3_tabn+0.1, he3_pabn, 'he3\_pabn,');
  text(he3_tabn+0.1, he3_pabn-1.5, 'he3\_tabn');
  text(1.3,10,   'he3\_tc(p)');
  text(1.7,26, 'he3\_tab(p)');

  text(2.17,30, 'A','fontweight','bold');
  text(1.3,20, 'B','fontweight','bold');
  text(2,5, 'Normal','fontweight','bold');
  text(2.1,40, 'Solid','fontweight','bold');
  xlabel('temperature, mK');
  ylabel('pressure, bar');

  print phase1.eps -deps -color

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clf; hold on;
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

  text(he3_tsmin, 15, 'he3\_psmin,');
  text(he3_tsmin,  7, 'he3\_tsmin');
  text(he3_tcr*1.5, 0.15, 'he3\_pcr,');
  text(he3_tcr*1.5, 0.07, 'he3\_tcr');
  text(0.005,80, 'he3\_pmelt(t)');
  text(0.6,1e-4,  'he3\_pvap(t)');

  text(0.01,0.01, 'Liquid','fontweight','bold');
  text(0.01,1000, 'Solid','fontweight','bold');
  text(5,0.003, 'Gas','fontweight','bold');
  text(0.0012,10, 'B','fontweight','bold');

  xlim([min(temp) max(temp)]);
  ylim([10^-5 10^4]);
  xlabel('temperature, K');
  ylabel('pressure, bar');

  print phase2.eps -deps -color




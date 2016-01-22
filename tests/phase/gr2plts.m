#!/usr/bin/octave -qf

  addpath ../../matlab

  figure; clf; hold on;
  t=0.5:0.01:3;
  p=34.32:0.001:34.4;

  % greywall and PLTS points:
  tg = [he3_ta he3_tb he3_ts];
  pg = [he3_pa he3_pb he3_ps];
  tp = [he3_ta_plts he3_tb_plts he3_ts_plts];
  pp = [he3_pa_plts he3_pb_plts he3_ps_plts];


  % conversion Greywall -> PLTS
  p1p1 = sum(pg.*pg);
  t1t1 = sum(tg.*tg);
  t1p1 = sum(tg.*pg);
  ddet = t1t1*p1p1 - t1p1*t1p1;
  t2t1 = sum(tp.*tg);
  t2p1 = sum(tp.*pg);
  p2t1 = sum(pp.*tg);
  p2p1 = sum(pp.*pg);

  % t2 = c(1)*t1 + c(2)*p1
  % p2 = c(3)*t1 + c(4)*p1
%  c(1) = (p1p1*t2t1 - t1p1*t2p1)/ddet;
%  c(2) = (t1t1*t2p1 - t1p1*t2t1)/ddet;
%  c(3) = (p1p1*p2t1 - t1p1*p2p1)/ddet;
%  c(4) = (t1t1*p2p1 - t1p1*p2t1)/ddet;
%  gr2plts_t = @(tg, pg) tg*c(1) + pg*c(2);
%  gr2plts_p = @(tg, pg) tg*c(3) + pg*c(4);

  % t2 = c(1)*t1 + c(2)*p1
  % p2 = c(3)*p1
%  c(1) = (p1p1*t2t1 - t1p1*t2p1)/ddet - 1;
%  c(2) = (t1t1*t2p1 - t1p1*t2t1)/ddet;
%  c(3) = p2p1/p1p1 - 1;
%  gr2plts_t = @(tg, pg) tg*(c(1)+1) + pg*c(2);
%  gr2plts_p = @(tg, pg) pg*(c(3)+1);

  % best method = 
  cp = polyfit(pg,pp./pg-1,0);
  ct = polyfit(tg,tp./tg-1,2);
  gr2plts_t = @(tg,pg) tg.*(polyval(ct,tg)+1);
  gr2plts_p = @(tg,pg) pg.*(polyval(cp,pg)+1);

  cp
  ct

  tabn = gr2plts_t(2.273,21.22)
  pabn = gr2plts_p(2.273,21.22)

%  c
  % Plot Greywall and PLTS
  plot(t, he3_pmelt(t/1000), 'b-');
  plot(t, he3_pmelt_plts(t/1000), 'r-');
  plot(tp, pp, 'ro');
  plot(tg, pg, 'bo');

  plot(he3_tab(p), p, 'b-');
  plot(he3_tc(p), p, 'b-');

  % plot converted:
  plot(gr2plts_t(tg,pg), gr2plts_p(tg,pg), 'g*');
  plot(gr2plts_t(t,he3_pmelt(t/1000)), gr2plts_p(t,he3_pmelt(t/1000)), 'g-');
  plot(gr2plts_t(he3_tab(p),p), gr2plts_p(he3_tab(p),p), 'g-');
  plot(gr2plts_t(he3_tc(p),p), gr2plts_p(he3_tc(p),p), 'g-');

  xlim([min(t) max(t)]);
  ylim([min(p) max(p)]);

  xlabel('terature, mK');
  ylabel('pure, bar');




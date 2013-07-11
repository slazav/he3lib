function plot_phase_diag()
  find_figure('He3 phase diagram'); clf;
  subplot(2,2,1); hold on;
  temp=0.3:0.01:3;
  plot(temp, he3_pmelt(temp), 'b-');
  press=0:0.1:36;
  plot(he3_tab(press), press, 'b-');
  plot(he3_tc(press), press, 'b-');
  plot(he3_tabn, he3_pabn, 'r*');
  plot(he3_tc(he3_pa), he3_pa, 'r*');
  xlabel('temperature, mK');
  ylabel('pressure, bar');

  subplot(2,2,2); hold on;
  temp=0.3:0.1:500;
  plot(temp, he3_pmelt(temp), 'b-');
  plot(he3_tc(press), press, 'b-');
  plot(he3_tab(press), press, 'b-');
  plot(he3_tabn, he3_pabn, 'r*');
  plot(he3_tc(he3_pa), he3_pa, 'r*');
  plot(he3_tsmin, he3_psmin, 'r*');
  plot(temp, he3_pvap(temp/1000), 'b-');
  xlabel('temperature, mK');
  ylabel('pressure, bar');

  subplot(2,2,3); hold on; title('Melting pressure');
  t=0.1:0.01:31;
  plot(t, he3_pmelt(t*1000), 'b-');
  plot(t, he3_pvap(t), 'b-');
  plot(he3_tsmin/1000, he3_psmin, 'r*');
  xlabel('temperature, K');
  ylabel('pressure, bar');

  subplot(2,2,4); hold on; title('Vapor pressure');
  temp=0.1:0.01:3.4;
  plot(temp, he3_pvap(temp), 'b-');
  xlabel('temperature, K');
  ylabel('pressure, bar');

end
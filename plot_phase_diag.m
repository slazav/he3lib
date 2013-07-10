function plot_phase_diag()
  addpath ./matlab
  find_figure('He3 phase diagram'); clf;
  subplot(2,1,1); hold on;
  temp=0.3:0.01:3;
  plot(temp, he3_pmelt(temp), 'b-');
  press=0:0.1:36;
  plot(he3_tc(press), press, 'r-');
%  plot(he3_tab(press), press, 'g-');
  plot(he3_tabn, he3_pabn, 'r*');
%  plot(he3_tc(he3_pa), he3_pa, 'r*');
  xlabel('temperature, mK');
  ylabel('pressure, bar');

  subplot(2,1,2); hold on;
  temp=0.3:0.1:250;
%  plot(temp, he3_pmelt(temp), 'b-');
%  plot(he3_tc(press), press, 'b-');
%  plot(he3_tab(press), press, 'b-');
  xlabel('temperature, mK');
  ylabel('pressure, bar');
end
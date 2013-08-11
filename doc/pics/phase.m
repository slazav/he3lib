#!/usr/bin/octave -qf


function title1(t)
  text(0.5,0.85, t, ...
    'units', 'normalized', 'horizontalalignment', 'center');
end

  addpath ../../matlab

  figure; clf;
  subplot(2,2,1); hold on;
  temp=0.3:0.01:3.7;
  plot(temp, he3_pmelt(temp/1000), 'b-');
  press=0:0.1:36;
  plot(he3_tab(press), press, 'b-');
  plot(he3_tc(press), press, 'b-');
  plot(he3_tabn, he3_pabn, 'ro');
  plot(he3_ta, he3_pa, 'ro');
  plot(he3_tb, he3_pb, 'ro');
  plot(he3_tneel, he3_pneel, 'go');
  xlim([0.5 3.7]);
  ylim([0 45]);
  text(he3_tneel-0.2, he3_pneel+6, 'he3\_pneel,', 'fontsize', 8);
  text(he3_tneel-0.2, he3_pneel+3, 'he3\_tneel', 'fontsize', 8);
  text(he3_tb-0.2, he3_pb+6, 'he3\_pb,', 'fontsize', 8);
  text(he3_tb-0.2, he3_pb+3, 'he3\_tb', 'fontsize', 8);
  text(he3_ta, he3_pa+6, 'he3\_pa,', 'fontsize', 8);
  text(he3_ta, he3_pa+3, 'he3\_ta', 'fontsize', 8);
  text(he3_tabn+0.1, he3_pabn, 'he3\_pabn,', 'fontsize', 8);
  text(he3_tabn+0.1, he3_pabn-3, 'he3\_tabn', 'fontsize', 8);
  text(1.0,10,   'he3\_tc(p)', 'fontsize', 8);
  text(1.2,26, 'he3\_tab(p)', 'fontsize', 8);

  text(2.15,29, 'A','fontweight','bold');
  text(1.3,18, 'B','fontweight','bold');
  text(2.7,27, 'Normal','fontweight','bold');
  text(3.1,37, 'Solid','fontweight','bold');
  xlabel('temperature, mK');
  ylabel('pressure, bar');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,2); hold on;
  temp=10.^(-3:0.1:2);
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

  text(he3_tsmin, 10, 'he3\_psmin,', 'fontsize', 8);
  text(he3_tsmin,  3, 'he3\_tsmin', 'fontsize', 8);
  text(he3_tcr*1.5, 0.3, 'he3\_pcr,', 'fontsize', 8);
  text(he3_tcr*1.5, 0.1, 'he3\_tcr', 'fontsize', 8);
  text(0.005,80, 'he3\_pmelt(t)', 'fontsize', 8);
  text(0.6,1e-4,  'he3\_pvap(t)', 'fontsize', 8);

  text(0.01,0.01, 'Liquid','fontweight','bold');
  text(0.01,1000, 'Solid','fontweight','bold');
  text(5,0.003, 'Gas','fontweight','bold');

  xlim([min(temp) max(temp)]);
  ylim([10^-5 10^4]);
  xlabel('temperature, K');
  ylabel('pressure, bar');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subplot(2,2,3); hold on; title1('Melting pressure');
  t=0.1:0.01:31;
  plot(t, he3_pmelt(t), 'b-');
  xlabel('temperature, K');
  ylabel('pressure, bar');

  subplot(2,2,4); hold on; title1('Vapor pressure');
  temp=0.1:0.01:3.4;
  plot(temp, he3_pvap(temp), 'b-');
  plot(he3_tcr, he3_pcr, 'ro');
  xlabel('temperature, K');
  ylabel('pressure, bar');

  print phase.eps -deps -color

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make table
f=fopen('phase.tex','w');
press=[0:3:33];
fprintf(f, '\\noindent\\begin{tabular}{r|');
for i = 1:length(press); fprintf(f, 'c'); end
fprintf(f, '}\n');

fprintf(f, 'P, bar')
for i = 1:length(press); fprintf(f, '& %2d', press(i)); end
fprintf(f, '\\\\\n');

fprintf(f, '\\hline');

fprintf(f, '$T_c$, mK')
for i = 1:length(press); fprintf(f, '& %5.3f', he3_tc(press(i))); end
fprintf(f, '\\\\\n');

fprintf(f, '$T_{AB}$, mK')
for i = 1:length(press);
  if (press(i)>he3_pabn) fprintf(f, '& %5.3f', he3_tab(press(i)));
  else  fprintf(f, '& --'); end
end
fprintf(f, '\\\\\n');

fprintf(f, '\\end{tabular}\n');
fclose(f);



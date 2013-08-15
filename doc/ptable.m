#!/usr/bin/octave -qf

function p1(f, x, t, fmt)
  fprintf(f, t)
  for i = 1:length(x);
    if isnan(x(i))
      fprintf(f, '& --');
    else
      fprintf(f, ['& %' fmt], x(i) );
    end
  end;
  fprintf(f, '\\\\\n');
end

addpath ../matlab

f=fopen('ptable.tex','w');
press=[0:3:33];
fprintf(f, '\\noindent\\begin{tabular}{r|');
for i = 1:length(press); fprintf(f, 'c'); end
fprintf(f, '}\n');

p1(f, press, 'P, bar', '2d');

fprintf(f, '\\hline');

p1(f, he3_tc(press),       '$T_c$, mK',      '5.3f');
tab=he3_tab(press); tab(find(press<he3_pabn)) = nan;
p1(f, tab,  '$T_{AB}$, mK',   '5.3f');

fprintf(f, '\\hline');

p1(f, he3_vm(press),       '$v_m$, cm$^3$/mol',      '5.2f');
p1(f, he3_c1(press)/100,   '$c_1$, m/s',             '5.1f');
p1(f, he3_gammaf(press),   '$\\gamma_f$, 1/(K mol)', '5.3f');
p1(f, he3_tmag(press)*1e3, '$T^\\star$, mK',         '5.1f');

fprintf(f, '\\hline\n');

p1(f, he3_rho(press),        '$\\rho$, g/cm$^3$',      '5.3f');
p1(f, he3_2n0(press)/1e38,   '$2N(0), 10^{38}$',       '5.3f');
p1(f, he3_pf(press)/1e-20,   '$p_F, 10^{-20}$ g cm/s', '5.3f');
p1(f, he3_meff(press)/1e-23, '$m^\\star, 10^{-23}$ g', '5.3f');
p1(f, he3_mm(press),         '$m^\\star/m_3$',         '5.3f');

p1(f, he3_vf(press)/1e2,     '$v_F$,~m/s',          '5.2f');
p1(f, he3_chi_n(press)/1e-9, '$\\chi_N, 10^{-9}$',  '5.1f');
p1(f, he3_f0s(press),        '$F_0^s$',             '5.2f');
p1(f, he3_f1s(press),        '$F_1^s$',             '5.2f');
p1(f, he3_f0a(press),        '$F_0^a$',             '5.2f');
p1(f, he3_f1a(press),        '$F_1^a$',             '5.2f');
p1(f, he3_a(press),          '$a$,~\\AA',           '5.3f');
p1(f, he3_gdk(press)/1e-9,   '$g_d/k_B$,~$\\mu$K',  '5.1f');
p1(f, he3_tfeff(press),      '$T_{F_{eff}}$,~K',    '5.3f');

fprintf(f, '\\hline\n');

p1(f, he3_tau_n0(1, press)/1e-6, '$\\tau_N(0,T_c)$, $\\mu$s', '5.3f');
p1(f, he3_scatt_l1a(press),    '$\\lambda_1^a$',      '5.3f');
p1(f, he3_scatt_g0(press),     '$\\gamma_0$',         '5.3f');
p1(f, he3_scatt_d0(press),     '$\\delta_0$',         '5.3f');
p1(f, he3_scatt_w0(press),     '$w_0$',               '5.3f');

fprintf(f, '\\hline\n');

p1(f, press, 'P, bar', '2d');

%%

fprintf(f, '\\end{tabular}\n');
fclose(f);



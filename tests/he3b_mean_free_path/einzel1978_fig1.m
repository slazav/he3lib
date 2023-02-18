#!/usr/bin/octave-cli -f
pkg load he3lib

t='Einzel JLTP32 (1978) p.35 fig.1'

% Einzel JLTP32 (1978) p.35 fig.1
addpath ../../matlab
figure; clf; hold on;
ttc = 0.00:0.01:1;
p=0;

gap=he3_trivgap(ttc,p);

tn0=he3_tau_n0(ttc,p);
tna=he3_tau_n_av(ttc,p);
t0=he3_tau0(ttc, p);
ta=he3_tau_av(ttc, p);

[x,y] = textread('einzel1978_fig1.dat', '%f %f', 'commentstyle', 'shell');
plot(x,y, 'k-');

plot(gap./ttc, tn0./t0, 'b-');
plot(gap./ttc, tna./ta, 'r-');
xlabel('\Delta/T')
ylabel('\tau_N/\tau')
legend('plot in the paper',...
       '\tau_N(0)/\tau(0)',...
       '<\tau_N(E)>/<\tau(E)>'...
)

title(t)
ylim([0 1]);
xlim([0 9]);

print -dpng -color einzel1978_fig1_.png


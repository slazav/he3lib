#!/usr/bin/octave -qf

addpath ../../matlab

figure; clf; hold on;

ttc=0:0.01:1.02;
p=0;

plot(ttc, he3_chi_b(ttc,00), 'r')
plot(ttc, he3_chi_b(ttc,10), 'g')
plot(ttc, he3_chi_b(ttc,20), 'b')
plot(ttc, he3_chi_b(ttc,30), 'm')
xlim([0 1.02]);
legend('0 bar', '10 bar', '20 bar', '30 bar');
text(0.2, 2e-8, '\chi_B', 'fontsize', 10);

xlabel('T/T_c');
ylabel('\chi_B/\chi_N');

print gap_sus.eps -deps "-S640,480" -color


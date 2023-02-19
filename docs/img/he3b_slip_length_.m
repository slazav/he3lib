#!/usr/bin/octave-cli -qf

pkg load he3lib

figure; clf; hold on;
xlabel('T/Tc');
ylabel('\zeta/l');

ttc=0:0.01:1.1;
v1=he3b_slip_length(ttc,0);
v2=he3b_slip_length(ttc,30);

plot(ttc,v1, 'r-');
plot(ttc,v2, 'b-');
legend('P=0 bar', 'P=30 bar');
xlim([0,1.1])

print he3b_slip_length.png -dpng "-S400,300"

#!/usr/bin/octave -qf

figure; hold on;
addpath ../../matlab

ttc=0.01:0.01:1;
gap=he3_bcsgap(ttc);
plot(ttc, he3_yosida(gap,ttc,0), 'r')
plot(ttc, he3_yosida(gap,ttc,2), 'g')
plot(ttc, he3_yosida(gap,ttc,4), 'b')

print yosida.eps -deps -color

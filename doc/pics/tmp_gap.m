#!/usr/bin/octave -qf

figure;
addpath ../../matlab

ttc=linspace(0,1,1e7);

g1=he3_bcsgap(ttc);
g2=he3_bcsgap_fast(ttc);

%plot(ttc, (g1-g2)./g1, 'r-');

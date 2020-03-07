#!/usr/bin/octave -qf

  addpath ../../matlab
  graphics_toolkit("gnuplot")
  figure; hold on;

function int = integrand(x,sgap,TTc,n)
    ksi = atanh(x)*TTc*2;
    Ek=sqrt(ksi.^2+sgap.^2);
    phi=(cosh(Ek/(2*TTc))).^(-2);
    int=(ksi/Ek).^n .* phi * 2*TTc .* cosh(ksi/2/TTc).^2;
end

ttc=[0.01 0.05 0.2 0.4 0.6 0.8 0.9 0.95 0.99];

x=0:0.001:0.999;

for i=1:length(ttc)
  t=ttc(i)
  gap=he3_bcsgap(t);
  plot(x, integrand(x,gap,t,0), 'r')
  plot(x, integrand(x,gap,t,2), 'g')
  plot(x, integrand(x,gap,t,4), 'b')
end


print yosida_int.eps -deps -color

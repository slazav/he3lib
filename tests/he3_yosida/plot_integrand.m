function plot_integrand()
% Plot yosida function integrand for different
% temperatures.
% It should not go to infinity at any T/Tc.
% This can be controlled by changing c in
% x = tanh(xi/c(ttc))
% c = 2 is a rather good value

  clf; hold on;

  x=0.001:0.001:0.999;

  for TTc=[0.05 0.5 0.99]
    gap=he3_bcsgap(TTc)
    I0 = yosida_int(x,gap,TTc,0);
    I1 = yosida_int(x,gap,TTc,1);
    I2 = yosida_int(x,gap,TTc,2);
    I3 = yosida_int(x,gap,TTc,3);
    I4 = yosida_int(x,gap,TTc,4);
    plot(x, I0/max(I0), 'r')
    plot(x, I1/max(I1), 'g')
    plot(x, I2/max(I2), 'b')
    plot(x, I3/max(I3), 'm')
    plot(x, I4/max(I4), 'c')
  end

%  print -deps -color plot_yosida_int.eps

end

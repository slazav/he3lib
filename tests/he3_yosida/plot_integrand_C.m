function plot_integrand()
% Plot yosida function integrand for different
% temperatures.
% It should not go to infinity at any T/Tc.
% This can be controlled by changing c in
% x = tanh(xi/c(ttc))
% c = 2 is a rather good value

  clf; hold on;

  x=0.001:0.001:0.999;

  for ttc=[0.05 0.5 0.99]
    gap=he3_bcsgap(ttc)
    dgap2=he3_bcsdgap2(ttc)
    I0 = integrand_c(x,gap,dgap2,ttc);
    plot(x, I0/max(I0), 'r')
  end

%  print -deps -color plot_integrand_c.eps

end

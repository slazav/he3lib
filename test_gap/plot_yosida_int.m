function plot_yosida_int()
% It should not go to infinity at any T/Tc.
% This can be controlled by changing c in
% x = tanh(xi/c(ttc))
% c = 2 is a rather good value

  function int = integrand(x,gap,TTc,n)
    c = 2;  % important power factor
    ksi = atanh(x) * c;
    Ek=sqrt(ksi.^2+gap.^2);
    phi=(cosh(Ek/(2*TTc))).^(-2) / 2 / TTc;
    int=(ksi/Ek).^n .* phi ./ (1-x.^2) * c;
  end

  figure;
  hold on;
  addpath ../matlab

  x=0.001:0.001:0.999;

  for TTc=[0.05 0.5 0.99]
    gap=he3_bcsgap(TTc)
    ksi = atanh(x)*TTc*2;
    Ek=sqrt(ksi.^2+gap.^2);
    I = integrand(x,gap,TTc,0);
    plot(x, I/I(1), 'r')
  end

  print -deps -color plot_yosida_int.eps

end

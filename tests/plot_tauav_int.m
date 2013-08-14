function plot_yosida_int()
% Plot integrand in tau_av and fpath calculation.
% It should not go to infinity at any T/Tc.
% This can be controlled by changing c in
% x = tanh(xi/c(ttc))
% c = 3 is a rather good value

  figure;
  hold on;
  addpath ../matlab

  p=0;
  g0 = he3_scatt_g0(p);
  d0 = he3_scatt_d0(p);

  function int = int_tauav(x,gap,ttc,n, g0, d0)
    c = 3;  % important power factor
    xi = atanh(x) * c;
    Ek=sqrt(xi.^2+gap.^2);
    phi=(cosh(Ek/(2*ttc))).^(-2);
    I = he3_coll_int(xi, ttc, gap, g0, d0);
    int = I .* phi ./ (1-x.^2) * c;
  end

  function int = int_fpath(x,gap,ttc,n, g0, d0)
    c = 2;  % important power factor
    xi = atanh(x) * c;
    Ek=sqrt(xi.^2+gap.^2);
    phi=(cosh(Ek/(2*ttc))).^(-2);
    fp = 1./(1 + exp(Ek./ttc));
    I = he3_coll_int(xi, ttc, gap, g0, d0);
    int = 1./I.^2 .* (xi./Ek).^2 .* fp ./ (1-x.^2) * c;
  end


  x=0.001:0.001:0.999;

  for TTc=[0.01 0.5 0.999]
    gap=he3_trivgap(TTc, p)
    ksi = atanh(x)*TTc*2;
    Ek=sqrt(ksi.^2+gap.^2);
    I = int_tauav(x,gap,TTc,0, g0, d0);
    plot(x, I/I(1), 'r')

    I = int_fpath(x,gap,TTc,0, g0, d0);
    plot(x, I/max(I), 'b')
  end

  print -deps -color plot_tauav_int.eps

end

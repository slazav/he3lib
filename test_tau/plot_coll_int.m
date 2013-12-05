function plot_coll_int()
% Plot collision integral

  figure;
  hold on;
  addpath ../matlab

  p=0;

  g0=he3_scatt_g0(p);
  d0=he3_scatt_d0(p);
  w0=1 - 2*g0/3 + d0;


  function int = integrand(x,gap,TTc,n)
    c = 2;  % important power factor
    ksi = atanh(x) * c;
    Ek=sqrt(ksi.^2+gap.^2);
    phi=(cosh(Ek/(2*TTc))).^(-2);
    int=(ksi/Ek).^n .* phi ./ (1-x.^2) * c;
  end


  ttc=0:0.01:1;
  for x=[0 0.8 0.99]
    gap=he3_trivgap(ttc,p);
    Y0=he3_yosida(ttc, gap, 0);

    xi = atanh(x);
    Ek=sqrt(xi.^2+gap.^2);

    I0 = 3/2/pi * gap./ttc .* he3_yosida(ttc,gap, 0D0) * w0;

    I = he3_coll_int(xi, ttc, gap, g0, d0);
    Il = he3_coll_int_lt(xi, ttc, gap, g0, d0);
    Ih = he3_coll_int_ht(xi, ttc, gap, g0, d0);

    semilogy(ttc, I, 'r')
    semilogy(ttc, Il, 'g')
    semilogy(ttc, Ih, 'b')
  end

  ylim(10.^[-2 1]);

  xlabel('T/T_c')
  ylabel('I')
  legend('full range approximation', 'low T', 'high T')

  print -deps -color plot_coll_int.eps

end

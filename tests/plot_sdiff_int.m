function plot_sdiff_int()
% Plot integrand for sdiff calculation.
% It should not go to infinity at any T/Tc.
% This can be controlled by changing c in
% x = tanh(xi/c(ttc))
% c = 3 is a rather good value

  figure;
  hold on;
  addpath ../matlab

  function int = int1(x, th, ttc, gap, o0, lambda, td, type)
    C = 3.5*ttc;  % power factor
    xi = atanh(x) * C;
    Ek=sqrt(xi.^2+gap.^2);
    phi=(cosh(Ek/(2*ttc))).^(-2) /2/ttc;
    u = xi./Ek;
    v = gap./Ek;
    kz=sin(th);
    kp=cos(th);
    o1 = o0*(1+lambda);

    Sm2 = 1 - kp.^2/2 .* (1-u.^2);
    Sp2 = u.^2 + (1-u.^2) .* kz.^2;
    t = td ./ (1 - i*o0*td);
    s = o1 .* t;

    if type==1
      int = t .* 0.5.*kp .* kp.^2 ...
        .* (Sm2 - i*Sp2.*s) ./ (1 + Sp2.*s.^2) ...
        .* phi .* C./(1-x.^2);
    elseif type==2
      int = t .* kp .* kz.^2 ...
        .* (Sm2 - i*Sp2.*s) ./ (1 + Sp2.*s.^2) ...
        .* phi .* C./(1-x.^2);
    elseif type==3
      int = td .* 0.5.*kp .* kp.^2 ...
        .* (Sm2 + u.^2 .* (o1.*td).^2) ./ (1 + Sp2.*(o1.*td).^2) ...
        .* phi .* C./(1-x.^2);
    elseif type==4
      int = td .* kp .* kz.^2 ...
        .* (Sm2 + u.^2 .* (o1.*td).^2) ./ (1 + Sp2.*(o1.*td).^2) ...
        .* phi .* C./(1-x.^2);
    end
  end


  p=0;
  x=0.01:0.01:0.99;

  for ttc=[0.01 0.05 0.999]
    gap=he3_trivgap(ttc, p)

    type=1;
    I = int1(x,89, gap, ttc, 600e3/2/pi, 1.1, he3_tau_dperp(ttc,p), type);
    plot(x, real(I)/real(I(1)), 'r')
    if type<3; plot(x, imag(I)/imag(I(1)), 'm'); end

    I = int1(x,45, gap, ttc, 600e3/2/pi, 1.1, he3_tau_dperp(ttc,p), type);
    plot(x, real(I)/real(I(1)), 'g')
    if type<3; plot(x, imag(I)/imag(I(1)), 'y'); end

    I = int1(x,0, gap, ttc, 600e3/2/pi, 1.1, he3_tau_dperp(ttc,p), type);
    plot(x, real(I)/real(I(1)), 'b')
    if type<3; plot(x, imag(I)/imag(I(1)), 'c'); end

  end

  print -deps -color plot_sdiff_int.eps

end

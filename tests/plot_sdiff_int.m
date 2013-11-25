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
    C = 3.5;  % power factor
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
        .* (Sm2 - i*Sp2.*s) ./ (1 + Sp2.*s.^2);
    elseif type==2
      int = t .* kp .* kz.^2 ...
        .* (Sm2 - i*Sp2.*s) ./ (1 + Sp2.*s.^2);
    elseif type==3
      int = td .* 0.5.*kp .* kp.^2 ...
        .* (Sm2 + u.^2 .* (o1.*td).^2) ./ (1 + Sp2.*(o1.*td).^2);
    elseif type==4
      int = td .* kp .* kz.^2 ...
        .* (Sm2 + u.^2 .* (o1.*td).^2) ./ (1 + Sp2.*(o1.*td).^2);
    elseif type==5
      A=1 + kz.^2.*s.^2;
      B=(1-kz.^2).*s.^2;
      int = t .* 0.5.*kp .* kp.^2 ...
        .* ((1 - kp.^2/2 .* (1-u.^2)) - i*(u.^2 + (1-u.^2) .* kz.^2).*s) ...
        ./A .* (1 - B./A.*u.^2);

%        .* ((1 - kp.^2/2 .* (1-u.^2)) - i*(u.^2 + (1-u.^2) .* kz.^2).*s) ...
%        ./ (1 + (u.^2 + (1-u.^2) .* kz.^2).*s.^2);
    end
    int = int .* phi .* C./(1-x.^2);
  end


  p=0;
  x=0.1:0.001:0.3;
  o0 = 600e1/2/pi;
    lambda=0.01;
  ttcs=[0.13 0.12];


  type = 1;
  th=0.1;
  for ttc=ttcs
    gap=he3_trivgap(ttc, p);
    I = int1(x, th, gap, ttc, o0, lambda, he3_tau_dperp(ttc,p), type);
    plot(x, real(I), 'b')
    plot(x, imag(I), 'r')
  end


  type = 5;
  for ttc=ttcs
    gap=he3_trivgap(ttc, p);
    I = int1(x, th, gap, ttc, o0, lambda, he3_tau_dperp(ttc,p), type);
    plot(x, real(I), 'c')
    plot(x, imag(I), 'm')

  end

  print -deps -color plot_sdiff_int.eps

end

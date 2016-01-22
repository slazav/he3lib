function plot_sdiff_int1()
% Plot integrand for sdiff calculation.
% New one, with 1D integration
% It should not go to infinity at any T/Tc.
% This can be controlled by changing c in
% x = tanh(xi/c(ttc))

  figure;
  hold on;
  addpath ../../matlab

  function int = int1(x, ttc, gap, o0, lambda, td, type)
    C = 4.0;  % power factor
    xi = atanh(x) * C;
    Ek=sqrt(xi.^2+gap.^2);
    phi=(cosh(Ek/(2*ttc))).^(-2) /2/ttc;
    uu = (xi./Ek).^2;
    vv = (gap./Ek).^2;

    o1 = o0*(1+lambda);
    t = td ./ (1 - i*o0*td);
    s = o1 .* t;

    if type == 1 || type == 2
      AA = (0.5 - i*s).*vv;
      BB = 0.5*(1+uu) - i*s.*uu;
      CC = s.^2 .* vv;
      DD = 1 + s.^2 .* uu;
    else
      AA = 0.5*vv;
      BB = 0.5*(1+uu) + (o1.*td).^2 .* uu;
      CC = (o1.*td).^2 .* vv;
      DD = 1 + (o1.*td).^2 .* uu;
    end

    AC = AA./CC;
    CD = CC./DD;
    D1 = BB./AA-DD./CC;


    T1 = sqrt(CD) .* atan(sqrt(CD));
    I1 = AC + (BB./CC - AA.*DD./CC.^2).*T1;
    I2 = AC/3.0 + AC.*D1.*(1-T1./CD);

    ii = find(CC<1e-4);
    I1(ii) = (AA(ii)/3 + BB(ii) - AA(ii).*CC(ii)./DD(ii)/5 - BB(ii).*CC(ii)./DD(ii)/3)./DD(ii);
    I2(ii) = (AA(ii)/5 + BB(ii)/3 - AA(ii).*CC(ii)./DD(ii)/7 - BB(ii).*CC(ii)./DD(ii)/5)./DD(ii);

    if type==1
      int = t .* (I1-I2)/2.0;
    elseif type==2
      int = t * I2;
    elseif type==3
      int = td .* (I1-I2)/2.0;
    elseif type==4
      int = td .* I2;
    end
    int = int .* phi .* C./(1-x.^2);
  end


  p=5.2;
  x=0.1:0.001:0.3;
  o0 = 1/2/pi;
  lambda=0.01;
  ttcs=[0.13];


  type = 1;
  x=0:0.01:1;
  for ttc=ttcs
    gap=he3_trivgap(ttc, p);
    I = int1(x,  gap, ttc, o0, lambda, he3_tau_dperp(ttc,p), type);
    plot(x, real(I)/max(real(I)), 'b')
    plot(x, imag(I)/max(imag(I)), 'r')
  end

  type = 2;
  for ttc=ttcs
    gap=he3_trivgap(ttc, p);
    I = int1(x, gap, ttc, o0, lambda, he3_tau_dperp(ttc,p), type);
    plot(x, real(I)/max(real(I)), 'c')
    plot(x, imag(I)/max(imag(I)), 'm')

  end

  print -deps -color plot_sdiff_int1.eps

end

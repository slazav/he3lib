function interp_y0()
  t='interpolation of Y0 function'

  addpath ../../matlab
  figure; clf; hold on;

  pp = fminunc(@minfunc, [0.39822 1])
%  pp = [0.39822 0.81178]

  ttc = 0.01:0.01:1;
  p = 30;
  gap=he3_trivgap(ttc,p);
  Y0  = he3_yosida(ttc,gap,0);
  Yi  = y0int(ttc,gap, pp);

  Ylt = sqrt(2*pi*gap./ttc) .* exp(-gap./ttc) .* (1+0.375*ttc./gap);

  plot(ttc, Y0./Ylt, 'r-');
  plot(ttc, Yi, 'g-');

  xlim([0 1]);
  ylim([0.4 1.1]);
  xlabel('ttc')
  ylabel('Y0/Y0lt')
  title(t);
  legend('Y0', 'Y0 interpolation')
  print -deps -color interp_y0.eps
end

function Yi=y0int(ttc,gap,pp)
%  Yi = (1 - pp(1)*(gap(1)-gap).^0.8 );
  Yi = (1 - pp(1)*(1-gap/gap(1)).^0.5 + pp(2)*(gap(1)-gap) );
end

function F = minfunc(pp)
  p=0;
  ttc = 0.01:0.01:0.99;
  gap=he3_trivgap(ttc,p);

  Y0 = he3_yosida(ttc, gap, 0);
  Ylt = sqrt(2*pi*gap./ttc) .* exp(-gap./ttc) .* (1+0.375*ttc./gap);
  Yi = y0int(ttc, gap, pp);
  F = sum((Y0./Ylt - Yi).^2);
end

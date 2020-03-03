# interpolated functions according with Eizel-2003
#  octave version + check of fortran functions

function interp1_y0()
  t='interpolation of Y0 function'

  addpath ../../octave
  figure; clf; hold on;

  P=0;
  ttc=0.05:0.01:1;

  gap = he3_gap(ttc,P);

  plot(ttc, he3_yosida(ttc, gap,0), 'c-');
  plot(ttc, y0int(ttc,P), 'r-');
  plot(ttc, ySint(ttc,P), 'b-');
  plot(ttc, yCint(ttc,P), 'g-');
  plot(ttc, he3_yosida_0(ttc,P), 'r.');
  plot(ttc, he3_yosida_s(ttc,P), 'b.');
  plot(ttc, he3_yosida_c(ttc,P), 'g.');

  xlim([0 1]);
  xlabel('ttc')
  ylabel('Y0')
  title(t);
  legend('Y0', 'Y0 interpolation', 'Ys interpolation', 'Yc interpolation')
  print -deps -color interp_y0.eps
end



function y=y0int(ttc,P)
  gap = he3_gap(0,P);
  y00 = sqrt(2*pi*gap./ttc).*exp(-gap./ttc);
  y0  = y00.*(1+3/8.0*ttc./gap);
  y0tc = sqrt(2*pi*gap).*exp(-gap).*(1+3/8.0./gap);
  k = (2.5-gap)/(1- y0tc);
  #  k = 2.388693; for BCS gap
  y = y0 .* (1-ttc.^k) + exp(gap-gap./ttc).*ttc.^(k-0.5);
end


function y=ySint(ttc,P)
  gap = he3_gap(0,P);
  dcbcn = he3_dcbcn(P);
  y00 = sqrt(2*pi*gap./ttc).*exp(-gap./ttc);
  y0  = 3/pi^2 * gap./ttc .* y00 .* (1 + 15/8.0*ttc./gap);
  y0tc  = 3/pi^2 * gap .* sqrt(2*pi*gap).*exp(-gap) .* (1 + 15/8.0./gap);

  k = (dcbcn + 1.5 - gap)/(1- y0tc);
  #  k = 2.150244; for BCS gap
  y = y0 .* (1-ttc.^k) + exp(gap-gap./ttc).*ttc.^(k-1.5);
end

function y=yCint(ttc,P)
  gap = he3_gap(0,P);
  #gap = he3_bcsgap(P);
  dcbcn = he3_dcbcn(P);
  y00 = sqrt(2*pi*gap./ttc).*exp(-gap./ttc);
  y0  = 3*(gap./ttc/pi).^2 .* y00 .* (1 + 11/8.0*ttc./gap);
  y0tc  = 3/pi^2 * gap.^2 .* sqrt(2*pi*gap).*exp(-gap) .* (1 + 11/8.0./gap)

  ytc = 1+dcbcn
  l2 = dcbcn*(1-31/144*1.036*dcbcn^2)

  k = (3*l2/ytc + 2.5 - gap)/(1- y0tc/ytc)
  #  k = 2.811729; for BCS gap
  c = (ytc-1)/y0tc - 1;
  y = y0 .* (1+c*ttc.^k) + exp(gap-gap./ttc).*ttc.^(k-2.5);

end

function plot_sdiff_per_lt2()

  addpath ../matlab

  figure; clf; hold on;

  p=4.1;
  f=[1];
  freqs=[f*0.1 f f*10 f*1e2 f*1e3 f*1e4 f*1e5 f*1e6];
  freqs=1e6;
  ttcr=5:0.1:10;
  ttc=1./ttcr;

  vf=he3_vf(p);
  f0a = he3_f0a(p);
  gap=he3_trivgap(ttc,p);
  y0=he3_yosida(ttc,gap,0);
  y2=he3_yosida(ttc,gap,2);
  y4=he3_yosida(ttc,gap,4);
%  chi = (2 + y0) / (3 + f0a*(2 + y0));
  chi = 2 / (3 + f0a*2);
  l = -f0a*chi;

  kxx = pi/8 * l*(2+l)^3/(1+l)^5;
  k0xx = kxx * 1
  k2xx = kxx * 1/2 * (l^2+2*l+4)
  k4xx = kxx * 3/8 * (8 + l*(2+l)*(4+2*l+l^2))

  kzz = pi/4 * (2+l)^2/(1+l)^5;
  k0zz = kzz * 1
  k2zz = -kzz * 1/2 * (l^2+2*l-2)
  k4zz = -kzz * 1/8 * (-8+l*(2+l)*(8+2*l+l^2))

  kk0xx=k0xx*2*gamma((0+1)/2)*2^(-0.5) / chi
  kk2xx=k2xx*2*gamma((2+1)/2)*2^(0.5)  / chi
  kk4xx=k4xx*2*gamma((4+1)/2)*2^(1.5)  / chi

  kk0zz=k0zz*2*gamma((0+1)/2)*2^(-0.5) / chi
  kk2zz=k2zz*2*gamma((2+1)/2)*2^(0.5)  / chi
  kk4zz=k4zz*2*gamma((4+1)/2)*2^(1.5)  / chi

  t0=(ttc./gap).^(-0.5);
  t2=(ttc./gap).^(0.5);
  t4=(ttc./gap).^(1.5);

  for f=freqs;

    DX = vf^2./(2*pi*f) .* exp(-gap./ttc);

    Dxx4= DX .* (kk0xx*t0 + kk2xx*t2 + kk4xx*t4);
    Dzz4= DX .* (kk0zz*t0 + kk2zz*t2 + kk4zz*t4);

    semilogy(ttcr, he3_diff_perp_xx(ttc, p, f), 'b-', 'linewidth', 2);
    semilogy(ttcr, he3_diff_perp_zz(ttc, p, f), 'r-', 'linewidth', 2);

    semilogy(ttcr, Dxx4, 'k-');
    semilogy(ttcr, Dzz4, 'k-');
  end
%  semilogy(ttcr, he3_diff_hperp_zz(ttc, p), 'g-');

  fixaxis;

  xlabel('T_c/T');
  ylabel('D_\perp');

  print -deps -color plot_sdiff_per_lt.eps
end
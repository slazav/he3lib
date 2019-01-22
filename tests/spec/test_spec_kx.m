function test_spec()
  addpath ../matlab
  figure(1); clf; hold on

  ttc=0;
  P=0;
  km = 5000;
  wm = 2*pi*2e6;

  fB = he3_nu_b(ttc,P);
  f0 = 500000;

  ak=0;
  bk=90/180*pi

  an=0/180*pi
  bn=30/180*pi


  H = 2*pi*f0/he3_gyro
  kv=linspace(0,km, 100);
  wv=linspace(1e-6, wm, 100);

  k=1/(2*pi*1e6);
  plot(kv, he3b_spec1(ttc,P,H, kv, ak,bk, an,bn)*k, 'r-', 'linewidth', 2)
  plot(kv, he3b_spec2(ttc,P,H, kv, ak,bk, an,bn)*k, 'r-', 'linewidth', 2)
  plot(kv, he3b_spec3(ttc,P,H, kv, ak,bk, an,bn)*k, 'r-', 'linewidth', 2)

  k2a = he3b_spec_kx2a(ttc,P,H, wv, an,bn);
  k2b = he3b_spec_kx2b(ttc,P,H, wv, an,bn);
  k2c = he3b_spec_kx2c(ttc,P,H, wv, an,bn);

  wa = wv;
  za = find(k2a(1:end-1)<0 & k2a(2:end)> 0);
  wa(za) = wa(za) - (wa(za+1)-wa(za))/(k2a(za+1)-k2a(za))*k2a(za);
  k2a(za) = 0;

  [wa, k2a] = fixk(wv, k2a);
  [wb, k2b] = fixk(wv, k2b);
  [wc, k2c] = fixk(wv, k2c);

  iia = find(k2a >=0 & k2a < km^2);
  iib = find(k2b >=0 & k2b < km^2);
  iic = find(k2c >=0 & k2c < km^2);

  plot(sqrt(k2a(iia)), wa(iia)*k, 'g-', 'linewidth', 2)
  plot(sqrt(k2b(iib)), wb(iib)*k, 'c-', 'linewidth', 2)
  plot(sqrt(k2c(iic)), wc(iic)*k, 'm-', 'linewidth', 2)

  plot([0 km], [f0 f0]/1e6, 'k-')
  plot([0 km], [fB fB]/1e6, 'b-')
  text(km, f0/1e6, 'f0')
  text(km, fB/1e6, 'fB')
  xlabel('k, 1/cm')
  ylabel('f, MHz')
end


function [w k2] = fixk(w, k2)
  za = find(k2(1:end-1)<=0 & k2(2:end)> 0);
  w(za) = w(za) - (w(za+1)-w(za))/(k2(za+1)-k2(za))*k2(za);
  k2(za) = 0;

  za = find(k2(1:end-1)>0 & k2(2:end)<=0);
  w(za+1) = w(za+1) - (w(za)-w(za+1))/(k2(za)-k2(za+1))*k2(za+1);
  k2(za+1) = 0;

end
function test_spec()
  figure(1); clf; hold on

  ttc=0.1;
  P=0;
  km = 2000;

  f0 = 300e3;
  fB = he3_nu_b(ttc,P);

  H = 2*pi*f0/he3_gyro;


  kv=linspace(0,km, 1000);
  ak=0;
  bk=pi/2;

  an=0/180*pi;
  bn=10/180*pi;

  k=1/(2*pi*1e6);
  plot(kv, he3b_spec1(ttc,P,H, kv, ak,bk, an,bn)*k, 'r-', 'linewidth', 2)
  plot(kv, he3b_spec2(ttc,P,H, kv, ak,bk, an,bn)*k, 'r-', 'linewidth', 2)
  plot(kv, he3b_spec3(ttc,P,H, kv, ak,bk, an,bn)*k, 'r-', 'linewidth', 2)

  plot(kv, he3b_spec1s(ttc,P,H, kv, ak,bk, an,bn)*k, 'b--')
  plot(kv, he3b_spec2s(ttc,P,H, kv, ak,bk, an,bn)*k, 'b--')
  plot(kv, he3b_spec3s(ttc,P,H, kv, ak,bk, an,bn)*k, 'b--')
  plot([0 km], [f0 f0]/1e6, 'k-')
  plot([0 km], [fB fB]/1e6, 'b-')
  text(km, f0/1e6, 'f0')
  text(km, fB/1e6, 'fB')
  xlabel('k, 1/cm')
  ylabel('f, MHz')
end
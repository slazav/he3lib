function test_spec()
  figure(1); clf; hold on

  ttc=0.1;
  P=0;
  H=255;

  kv=linspace(0,1e4,1000);
  ak=0;
  bk=10;

  an=20;
  bn=10;

  plot(kv, he3b_spec1(ttc,P,H, kv, ak,bk, an,bn), 'r-', 'linewidth', 2)
  plot(kv, he3b_spec2(ttc,P,H, kv, ak,bk, an,bn), 'r-', 'linewidth', 2)
  plot(kv, he3b_spec3(ttc,P,H, kv, ak,bk, an,bn), 'r-', 'linewidth', 2)

  plot(kv, he3b_spec1s(ttc,P,H, kv, ak,bk, an,bn), 'b--')
  plot(kv, he3b_spec2s(ttc,P,H, kv, ak,bk, an,bn), 'b--')
  plot(kv, he3b_spec3s(ttc,P,H, kv, ak,bk, an,bn), 'b--')
end
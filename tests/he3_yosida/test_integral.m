function test_integral()
  % integrate Yosida functions and compare result with he3lib:

  hold on;
  addpath ../../octave

  x=0.001:0.001:0.999;

  ttc=0.01:0.01:1;
  gap=he3_bcsgap(ttc);
  for n = 1:length(ttc)
    ff = @(x) yosida_int(x,gap(n),ttc(n),0);
    [y0(n), ier, nfun, err] = quad (ff, 0, 1);
  end

  plot(ttc, y0, 'r-')
  plot(ttc, he3_yosida(ttc, gap,0), 'b.')


end


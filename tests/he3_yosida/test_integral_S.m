function test_integral_S()
  % integrate Yosida_S function and compare result with he3lib:

  clf; hold on;
  addpath ../../octave

  ttc=0.01:0.01:1;
  for P = [0,10,20,30]
    gap=he3_trivgap(ttc,P);
    for n = 1:length(ttc)
      ff = @(x) integrand_s(x,gap(n),ttc(n));
      [yc(n), ier, nfun, err] = quad (ff, 0, 1);
    end
    yc = yc*3D0/const_pi**2;
    plot(ttc, yc, 'r-')
    plot(ttc, he3_yosida_s(ttc, gap), 'b.')
  end

  legend('he3lib', 'integration', ...
    'location', 'northwest')

end


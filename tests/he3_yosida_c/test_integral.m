function test_integral()
  % integrate Yosida_C function and compare result with he3lib:

  clf; hold on;
  addpath ../../octave

  for P = [0,10,20,30]
    ttc=0.01:0.01:1;
    gap=he3_trivgap(ttc,P);
    dgap2=he3_trivdgap2(ttc,P);
    for n = 1:length(ttc)
      ff = @(x) yosida_c_int(x,gap(n),dgap2(n),ttc(n));
      [yc(n), ier, nfun, err] = quad (ff, 0, 1);
    end
    yc = yc*3D0/const_pi**2
    plot(ttc, yc, 'r-')
    plot(ttc, he3_yosida_c(ttc, gap, dgap2), 'b.')
  end


end


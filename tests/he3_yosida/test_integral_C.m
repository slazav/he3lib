function test_integral_C()
  % integrate Yosida_C function and compare result with he3lib:

  clf; hold on;
  addpath ../../octave

  for P = [0,10,20,30]
    ttc=0.01:0.01:1;
    gap=he3_trivgap(ttc,P);
    dgap2=he3_trivdgap2(ttc,P);
    for n = 1:length(ttc)
      ff = @(x) integrand_c(x,gap(n),dgap2(n),ttc(n));
      [yc(n), ier, nfun, err] = quad (ff, 0, 1);
    end
    yc = yc*3D0/const_pi**2;
    plot(ttc, yc, 'r-')
    plot(ttc, he3_yosida_c(ttc, gap, dgap2), 'b.')

    % from entropy derivative
    ttc=0.01:0.01:0.99;
    dt = 0.01;
    yCs = ...
       he3_yosida_s(ttc, he3_trivgap(ttc,P)) + ttc.*( ...
       he3_yosida_s(ttc+dt, he3_trivgap(ttc+dt,P)) - 
       he3_yosida_s(ttc-dt, he3_trivgap(ttc-dt,P)))/2/dt;
    plot(ttc, yCs, 'g.')

  legend('he3lib', 'integration','entropy derivative', ...
    'location', 'northwest')

  end


end


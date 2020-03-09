function test_heat_cap()
  % integrate Yosida functions and compare result with he3lib:

  clf; hold on;
  addpath ../../octave

  for P = [0,10,20,30]
    ttc=0.01:0.01:1;
    gap=he3_trivgap(ttc,P);
    dgap2=he3_trivdgap2(ttc,P);
    plot(ttc, he3_yosida_c(ttc, gap, dgap2)./ttc.^2, 'b-')
  end


end


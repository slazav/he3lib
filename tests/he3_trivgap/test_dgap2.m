function test_dgap2()

  clf; hold on;
  addpath ../../octave

  dt  = 0.01;
  ttc = 0.01:0.01:0.99;

  for p=[0,10,20,30]
    gap = he3_trivgap(ttc, p);

    dgap2 = (he3_trivgap(ttc+dt,p).^2-he3_trivgap(ttc-dt,p).^2)/(2*dt);

    dgap2l = he3_trivdgap2(ttc,p);

    plot(ttc, dgap2, 'r-')
    plot(ttc, dgap2l, 'b.')

  end

  plot(ttc, he3_bcsdgap2(ttc), 'b-')

end


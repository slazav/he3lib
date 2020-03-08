function test_dgap2()
  % check he3_bcsgap2

  clf; hold on;
  addpath ../../octave

  dt  = 0.01;
  ttc = 0.01:0.01:0.99;
  gap = he3_bcsgap(ttc);

  dgap2 = (he3_bcsgap(ttc+dt).^2-he3_bcsgap(ttc-dt).^2)/(2*dt);

  dgap2l = he3_bcsdgap2(ttc);

  plot(ttc, dgap2, 'r-')
  plot(ttc, dgap2l, 'b.')


end


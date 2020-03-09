function test_gap()

  clf; hold on;
  addpath ../../octave

  dt  = 0.01;
  ttc = 0.01:0.01:0.99;

  for p=[0,10,20,30]
    gap = he3_trivgap(ttc, p);
    plot(ttc, gap, 'r-')
  end

  plot(ttc, he3_bcsgap(ttc), 'b-')


end


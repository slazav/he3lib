function plot_vm()

  t = 0:0.01:4.1;
  vm = he4_pmelt(t)
  figure; clf; hold on;
  plot(t, vm, 'r-')

end
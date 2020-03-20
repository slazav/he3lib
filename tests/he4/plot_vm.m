function plot_vm()

  t = 2.1:0.001:2.3;
  vm = he4_vm(t)
  figure; clf; hold on;
  plot(t, vm, 'r-')

end
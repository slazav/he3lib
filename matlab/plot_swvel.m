function plot_flegg()
  find_figure('Spin wave velocity'); clf;
  subplot(2,1,1); hold on;
  for pres=0:5:36;
    temp=0:0.01:1;
    plot(temp, he3_swvel_par(pres, temp), 'r-');
  end
  xlabel('temperature, mK');
  ylabel('Spin wave velocity, cm/s');

  subplot(2,1,2); hold on;
  pres=0:1:36;
  plot(pres, he3_swvel(pres, 0), 'r-');
  plot(pres, he3_swvel_par(pres, 0), 'b-');
  plot(pres, he3_swvel_per(pres, 0), 'g-');
  legend('C', 'Cpar', 'Cperp');
end
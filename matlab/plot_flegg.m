function plot_flegg()
  find_figure('Leggett freq'); clf; hold on;
  for pres=0:1:18;
    temp=0:0.01:1;
    plot(temp, he3_flegg(pres, temp), 'r-');
  end
  for pres=19:1:36;
    temp=0.25:0.01:1;
    plot(temp, he3_flegg(pres, temp), 'b-');
  end
  xlabel('temperature, mK');
  ylabel('Leggett freq^2, Hz^2');
end
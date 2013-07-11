function plot_diff()
  find_figure('Spin diffusion'); clf; hold on;
  title('Spin diffusion');

  pres=0:0.1:36;
  temp=0:0.01:5;
  pp = repmat(pres', 1, length(temp));
  tt = repmat(temp, length(pres), 1);
  sigproc.plot_3di(temp,pres, he3_d_exp(pp, tt)');
  plot(temp, he3_pmelt(temp), 'b-');
  plot(he3_tab(pres), pres, 'b-');
  plot(he3_tc(pres), pres, 'b-');

  xlabel('temperature, mK');
  ylabel('pressure, bar');
end
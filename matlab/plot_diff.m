function plot_diff()
  find_figure('Spin diffusion'); clf; hold on;
  for pres=0:1:20;
    temp=0:0.01:5;
    plot(temp, he3_d_exp(pres, temp), 'g-');
%    plot(temp, he3_ds_exp(pres, temp), 'r-');
%    plot(temp, he3_dn_exp(pres, temp), 'b-');
  end
  xlabel('temperature, mK');
  ylabel('Leggett freq^2, Hz^2');
end
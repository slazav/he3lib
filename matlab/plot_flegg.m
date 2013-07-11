function plot_flegg()
  find_figure('Leggett freq'); clf; hold on;
  title('Leggett freq^2, Hz^2');
%  for pres=0:1:18;
%    temp=0:0.01:1;
%    plot(temp, he3_flegg(pres, temp), 'r-');
%  end
%  for pres=19:1:36;
%    temp=0.25:0.01:1;
%    plot(temp, he3_flegg(pres, temp), 'b-');
%  end
%  xlabel('temperature, mK');
%  ylabel('Leggett freq^2, Hz^2');

  pres=0:0.1:36;
  temp=0:0.01:5;
  pp = repmat(pres', 1, length(temp));
  tt = repmat(temp, length(pres), 1);
  sigproc.plot_3di(temp,pres, he3_flegg(pp, tt./he3_tc(pp))');
  plot(temp, he3_pmelt(temp), 'b-');
  plot(he3_tab(pres), pres, 'b-');
  plot(he3_tc(pres), pres, 'b-');
  xlabel('temperature, mK');
  ylabel('pressure, bar');
end
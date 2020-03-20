function plot_f()

  find_figure('transp'); clf; hold on;

  l=0:0.1:2;
  plot(l, f_eta(l), "b-")
  plot(l, f_k(l),  "r-")


end




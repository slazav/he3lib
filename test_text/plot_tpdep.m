function plot_tpdep(f)
  figure; clf; hold on;

  ttc = 0:0.01:1;
  p   = 5:5:25;

  plot(ttc, f(ttc,0), 'k.-');
  plot(ttc, f(ttc,30), 'r.-');

  N = length(p);
  for i=1:length(p)
    l=plot(ttc, f(ttc, p(i)), 'b.-');
  end


  legend('30 bar', '0 bar')
  xlabel('T/T_c');
end
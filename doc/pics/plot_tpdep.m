function plot_tpdep(f, lpos)

  ttc = 0:0.01:1;

  colors=[
   0.8 0.0 0.8
   1.0 0.0 0.0
   0.6 0.6 0.0
   0.0 0.6 0.0
   0.0 0.6 0.6
   0.0 0.0 1.0
   0.0 0.0 0.0
  ];
  p=0:5:30;

  for i=1:length(p)
    plot(ttc, f(ttc,  p(i)),...
       '-', 'color', colors(i,:), 'linewidth', 4);
  end

  if nargin>1
    legend('0 bar', '5 bar', '10 bar', '15 bar',...
           '20 bar', '25 bar', '30 bar', 'Location', lpos)
  end
end
function solve_cubic_test()

  figure(1); clf; hold on;

  p=[3 2 -3 0];
  p=[0 2 -3 1];
  p=[1 2 3 2];
  p=[1 2 -3 -1];

  [n, x1,x2,x3] = solve_cubic(p(1),p(2),p(3),p(4))

  xmin=min([x1 x2 x3 -1]);
  xmax=max([x1 x2 x3 1]);
  xw=xmax-xmin;

  x=linspace(xmin-xw/4,xmax+xw/4,100);

  plot(x, polyval(p,x), 'b-');

  plot(x1,0, 'r*')
  plot(x2,0, 'r*')
  plot(x3,0, 'r*')
end




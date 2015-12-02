function a()
  find_figure('trajectory'); clf; hold on
  [x,y,z,dx,dy,dz] = textread('a', '%f %f %f %f %f %f');
  k = 10;
  for i=1:length(x)
    x0=x(i);
    y0=y(i);
    z0=z(i);
    x1=x(i) + k*dx(i);
    y1=y(i) + k*dy(i);
    z1=z(i) + k*dz(i);
    plot([z0,z1], [y0,y1], 'r-')
  end
  plot(z,y, 'b.')
end
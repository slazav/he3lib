function plot_h()
% gap distortion at the AB boundary

  figure(1); clf; hold on;
  p=0:0.01:34;
  for H=0:200:6000
    t = he3_b2tab(p, H);
    ii=find(isfinite(t));
    t=t(ii); p=p(ii);
    t=[0 t 0]; p=[p(1) p p(end)];
    plot(t,p, 'r-');
  end

end

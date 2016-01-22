function plot_hab()
% gap distortion at the AB boundary

  figure(1); clf; hold on;
  for H=0:200:6000
    p=0:0.01:35;
    t = he3_b2tab(p, H).*he3_tc(p);
    ii=find(isfinite(t));
    if length(ii)==0; continue; end
    t=t(ii); p=p(ii);
    t=[0 t 0]; p=[p(1) p p(end)];
    plot(t,p, 'r-');
  end

  p=0:0.01:35;
  plot(he3_tc(p),p, 'b-');

  t=0:0.01:3;
  plot(t, he3_pmelt(t*1e-3), 'b-');
end

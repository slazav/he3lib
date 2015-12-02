function plot_b2gap()
% effective field

  figure(1); clf; hold on;
  p=0;
  ttc=-0.1:0.1:1.1;
  H=100;
  He = he3_b2heff(ttc,p,H);
  plot(ttc, He/H,'k')
  plot(ttc, ones(size(ttc))/(1+he3_f0a(p)),'k')
end

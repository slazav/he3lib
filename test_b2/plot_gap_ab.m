function plot_gap_ab()
% gap distortion at the AB boundary

  figure(1); clf; hold on;
  p=0;
  ttc=-0.1:0.01:1.1;
  Hc = he3_b2hcr(ttc,p);
  g0 = he3_gap(ttc,p);
  g1 = he3_b2gap1(ttc,p,Hc);
  g2 = he3_b2gap2(ttc,p,Hc);
  plot(ttc, g0,'k')
  plot(ttc, g1,'r')
  plot(ttc, g2,'b')
  title('Gap distortion at the AB boundary')
  legend('\Delta_0','\Delta_{perp}','\Delta_{par}')
end

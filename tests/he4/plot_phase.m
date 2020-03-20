function plot_phase()

  figure; clf; hold on;
  temp=0.01:0.01:6;

  semilogy(temp, he4_pmelt(temp), 'b-');
  semilogy(he4_tcr, he4_pcr, 'r*');
  semilogy(he4_tcv, he4_pcv, 'ro');
  semilogy(he4_tcm, he4_pcm, 'mo');

  press=linspace(he4_pcv, he4_pcm, 30);
  semilogy(he4_tc(press), press, 'r-');
  semilogy(temp, he4_pvap(temp), 'm-');

  xlabel('temperature, K');
  ylabel('pressure, bar');

  ylim([0.001 140])
#  print phase1.eps -deps "-S640,480" -color

end


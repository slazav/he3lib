#!/usr/bin/octave-cli -qf

  pkg load he3lib
  figure; clf;
  subplot(1,2,1); hold on;
  title("He3-He4 phase diagram")
  ylabel('temperature, K');
  xlabel('concentration');

  t=0:0.005:1;
  x=0:0.005:1;
  plot(he34_xdil(t), t, 'b-');
  plot(he34_xcon(t), t, 'b-');
  plot(x, he34_tlambda(x), 'b-');
  plot(he34_xcr, he34_tcr, 'ro');

  text(he34_xcr+0.05, 1.0, 'he34\_xcr');
  text(he34_xcr+0.05, 0.85, 'he34\_tcr');

  text(0.12, 0.1, 'he34\_xdil(T)');
  text(0.55, 0.3, 'he34\_xcon(T)');
  text(0.45, 1.5, 'he34\_tlambda(x)');

  subplot(1,2,2); hold on;
  title("pressure dependance")
  ylabel('temperature, K');
  xlabel('concentration');


  p=[0, 1, 2, 5, 10];
  c='bgcmk';

  for i=1:length(p)
    plot(he34_xdil_p(t,p(i)), t, [c(i), '-']);
    plot(he34_xcon_p(t,p(i)), t, [c(i), '-'], 'HandleVisibility','off');
    plot(x, he34_tlambda_p(x,p(i)),  [c(i), '-'], 'HandleVisibility','off');
  end
  xlim([0, 1.0])
  ylim([0, 1.4])

  p=0:1:30
  plot(he34_xcr_p(p), he34_tcr_p(p), 'r-', 'linewidth', 2);

  legend('0bar', '1bar', '2bar', '5bar', '10bar')


  print he34_phase.png -dpng "-S800,400"



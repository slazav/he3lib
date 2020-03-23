function test_he3he4

  x = 0:0.01:1;
  p = 0:2:10
  find_figure('test_he3he4'); clf; hold on;

  t = linspace(0, he34_tcr, 50);
  plot(x, he34_tlambda(x), 'r.-')
  plot(he34_xdil(t), t, 'r.-')
  plot(he34_xcon(t), t, 'r.-')
  plot(he34_xcr_p(p), he34_tcr_p(p) , 'b*-')

  for i=1:length(p)
    t = linspace(0, he34_tcr_p(p(i)), 50);

    plot(x, he34_tlambda_p(x,p(i)), 'b-')
    plot(he34_xdil_p(t,p(i)), t, 'b-')
    plot(he34_xcon_p(t,p(i)), t, 'b-')
  end


end
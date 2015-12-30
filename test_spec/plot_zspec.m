function plot1()
  % plot B-phase spin-wave spectra
  % for k=0 (uniform NMR)

  find_figure('CW NMR SPEC'); clf; hold on;

  for cb=cosd([10 45 85])
    %cb=0.99;           % \cos\beta_n
    a = 0.334:0.01:20; % omegaB/omegaL

    A = (1+a.^2)/2;
    B = a.^2*cb^2
    w1 = sqrt( A + sqrt(A.^2 - B));
    w2 = sqrt( A - sqrt(A.^2 - B));

    xcrd=1./a; % omegaL/OmegaB
    %xcrd=a;
    plot(xcrd,w1,'r-')
    plot(xcrd,w2,'b-')

    plot(xcrd, a*cb,'k--')
    plot(xcrd, 1 + a.^2/2 * (1-cb^2),'k--')
    plot(xcrd, ones(size(a)),'k-')
  end

  xlabel('omegaL/OmegaB')
  ylabel('omega/omegaL')
  ylim([0 2])

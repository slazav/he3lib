#!/usr/bin/octave -qf

 % Einzel-2003 interpolation of yosida function (BCS)
 % Good for understanding low- and high-temperature behaviour.

 % These are final formulas in the paper
 % For some reason they are not fully consistent with
 % the interpolation method (or I do not fully understood it).

  addpath ../../octave
  figure; clf; hold on;
  P = 30;
  ttc = 0:0.005:1;
  gap = he3_trivgap(ttc,P);
  dgap2 = he3_trivdgap2(ttc,P);
  y0 = he3_yosida(ttc, gap, 0);
  yS = he3_yosida_s(ttc, gap);
  yC = he3_yosida_c(ttc, gap, dgap2);


  function y = y0appr(ttc, P)
    gap0 = he3_trivgap(0,P); % WCP gap at T=0
    dcbcn = he3_dcbcn(P);
    dcbcn0 = 12D0/7D0/const_z3; % BCS

    % lowest-temperature approximation
    y00=sqrt(2*pi*gap0./ttc).*exp(-gap0./ttc);

    % next term in ttc/gap
    bet = 3.0/8.0;
    y0=y00.*(1 + bet*ttc./gap0);

    % value of low-temperature function at Tc:
    y00tc=sqrt(2*pi*gap0).*exp(-gap0);
    y0tc=y00tc.*(1 + bet/gap0);

    % slope near Tc
    s = 2 * dcbcn/dcbcn0;

    k= (s + 0.5 - gap0)/(1 - y0tc);
    y = y0.*(1-ttc.^k) + exp(gap0-gap0./ttc) .* ttc.^(k-0.5);
  end

  function y = ySappr(ttc, P)
    gap0 = he3_trivgap(0,P); % WCP gap at T=0
    dcbcn = he3_dcbcn(P);

    % lowest-temperature approximation
    y00 = 3/pi^2 * gap0./ttc .* sqrt(2*pi*gap0./ttc).*exp(-gap0./ttc);

    % next term in ttc/gap
    bet = 15.0/8.0;
    y0=y00.*(1 + bet*ttc./gap0);

    % value of low-temperature function at Tc:
    y00tc=3/pi^2 * gap0 * sqrt(2*pi*gap0).*exp(-gap0);
    y0tc=y00tc.*(1 + bet/gap0);

    % slope near Tc
    s = dcbcn;

    k= (s + 1.5 - gap0)/(1 - y0tc);
    y = y0.*(1-ttc.^k) + exp(gap0-gap0./ttc) .* ttc.^(k-1.5);
  end

  function y = yCappr(ttc, P)
    gap0 = he3_trivgap(0,P); % WCP gap at T=0
    dcbcn = he3_dcbcn(P);
    dcbcn0 = 12D0/7D0/const_z3; % BCS

    dt = 0.001;
    l2 = -(he3_trivdgap2(1,P)-he3_trivdgap2(1-dt,P))/dt/pi^2; % gap second derivative at Tc

    % lowest-temperature approximation
    y00 = 3*(gap0./ttc/pi).^2 .* sqrt(2*pi*gap0./ttc).*exp(-gap0./ttc);

    % next term in ttc/gap
    bet = 11.0/8.0;
    y0=y00.*(1 + bet*ttc./gap0);

    % value of low-temperature function at Tc:
    y00tc=3*(gap0/pi)^2 * sqrt(2*pi*gap0).*exp(-gap0);
    y0tc=y00tc.*(1 + bet/gap0);

    % value in Tc
    ytc = 1+dcbcn;

    % slope near Tc
    s = 3*l2 / ytc;


    k= (s/ytc + 2.5*gap0)/(gap0-y0tc/ytc);
    y = y0.*(1-ttc.^k) + ytc*exp(gap0-gap0./ttc) .* ttc.^(k-2.5);

  end

  y0a = y0appr(ttc, P);
  ySa = ySappr(ttc, P);
  yCa = yCappr(ttc, P);

  plot(ttc, y0, 'k-');
  plot(ttc, y0a, 'r-');

  plot(ttc, yS, 'k-');
  plot(ttc, ySa, 'b-');

  plot(ttc, yC, 'k-');
  plot(ttc, yCa, 'g-');

  plot(ttc, y0a./y0, 'r-');
  plot(ttc, ySa./yS, 'b-');
  plot(ttc, yCa./yC, 'g-');
  plot([0 1], [1 1], 'k-');

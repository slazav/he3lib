#!/usr/bin/octave -qf

 % Einzel-2003 interpolation of yosida functions (BCS)
 % Good for understanding low- and high-temperature behaviour.

 % These are final formulas in the paper
 % For some reason they are not fully consistent with
 % the interpolation method (or I do not fully understood it).
 % missing term in Yc formula is added.

  addpath ../../octave
  figure; clf; hold on;
  ttc = 0:0.005:1;
  gap = he3_bcsgap(ttc);
  dgap2 = he3_bcsdgap2(ttc);
  y0 = he3_yosida(ttc, gap, 0);
  yS = he3_yosida_s(ttc, gap);
  yC = he3_yosida_c(ttc, gap, dgap2);


  function y = y0appr(ttc)
    gap0 = pi/exp(const_euler); % BCS gap at T=0
    dcbcn0 = 12D0/7D0/const_z3;

    % lowest-temperature approximation
    y00=sqrt(2*pi*gap0./ttc).*exp(-gap0./ttc);

    % next term in ttc/gap
    bet = 3.0/8.0;
    y0=y00.*(1 + bet*ttc./gap0);

    % value of low-temperature function at Tc:
    y00tc=sqrt(2*pi*gap0).*exp(-gap0);
    y0tc=y00tc.*(1 + bet/gap0);

    % slope near Tc
    s = 2;

    k= (s + 0.5 - gap0)/(1 - y0tc);
    y = y0.*(1-ttc.^k) + exp(gap0-gap0./ttc) .* ttc.^(k-0.5);
  end

  function y = ySappr(ttc)
    gap0 = pi/exp(const_euler); % BCS gap at T=0
    dcbcn0 = 12D0/7D0/const_z3;

    % lowest-temperature approximation
    y00 = 3/pi^2 * gap0./ttc .* sqrt(2*pi*gap0./ttc).*exp(-gap0./ttc);

    % next term in ttc/gap
    bet = 15.0/8.0;
    y0=y00.*(1 + bet*ttc./gap0);

    % value of low-temperature function at Tc:
    y00tc=3/pi^2 * gap0 * sqrt(2*pi*gap0).*exp(-gap0);
    y0tc=y00tc.*(1 + bet/gap0);

    % slope near Tc
    s = dcbcn0;

    k= (s + 1.5 - gap0)/(1 - y0tc);
    y = y0.*(1-ttc.^k) + exp(gap0-gap0./ttc) .* ttc.^(k-1.5);
  end

  function y = yCappr(ttc)
    gap0 = pi/exp(const_euler); % BCS gap at T=0
    dcbcn0 = 12D0/7D0/const_z3;
    l2 = dcbcn0 * (1 - 31/144.0 * const_z5 * dcbcn0^2);

    % lowest-temperature approximation
    y00 = 3*(gap0./ttc/pi).^2 .* sqrt(2*pi*gap0./ttc).*exp(-gap0./ttc);

    % next term in ttc/gap
    bet = 11.0/8.0;
    y0=y00.*(1 + bet*ttc./gap0);

    % value of low-temperature function at Tc:
    y00tc=3*(gap0/pi)^2 * sqrt(2*pi*gap0).*exp(-gap0);
    y0tc=y00tc.*(1 + bet/gap0);

    % slope near Tc
    s = 3*l2 / (1+dcbcn0);

    % value in Tc
    ytc = 1+dcbcn0;

    k= (s + 2.5 - gap0)/(1 - y0tc/ytc);
    y = y0.*(1-ttc.^k) +ytc*exp(gap0-gap0./ttc) .* ttc.^(k-2.5);

  end


  y0a = y0appr(ttc);
  ySa = ySappr(ttc);
  yCa = yCappr(ttc);

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

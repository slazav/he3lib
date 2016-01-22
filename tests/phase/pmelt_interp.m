#!/usr/bin/octave -qf

  % interpolation of the melting pressure using various datasets

  addpath ../../matlab
  figure; clf; hold on;

  DGR = 10; % shift Greywall data

  % first define all available data:

  % greywall-1986 T = 0.0009 - 0.25 K
  function P = pmelt_gr(T)
   P = 34.3380e0 ...
     - 0.19652970e-1 * (T*1e3).^(-3) ...
     + 0.61880268e-1 * (T*1e3).^(-2) ...
     - 0.78803055e-1 * (T*1e3).^(-1) ...
     + 0.13050600e0 ...
     - 0.43519381e-1 * (T*1e3) ...
     + 0.13752791e-3 * (T*1e3).^2 ...
     - 0.17180436e-6 * (T*1e3).^3 ...
     - 0.22093906e-9 * (T*1e3).^4 ...
     + 0.85450245e-12* (T*1e3).^5;
  end

  % PLTS-2000 T = 0.0009 - 1 K
  function P = pmelt_plts(T)
    P =   (- 1.3855442e-12 * T.^(-3) ...
           + 4.5557026e-9  * T.^(-2) ...
           - 6.4430869e-6  * T.^(-1) ...
           + 3.4467434e0 ...
           - 4.4176438e0 * T.^1 ...
           + 1.5417437e1 * T.^2 ...
           - 3.5789858e1 * T.^3 ...
           + 7.1499125e1 * T.^4 ...
           - 1.0414379e2 * T.^5 ...
           + 1.0518538e2 * T.^6 ...
           - 6.9443767e1 * T.^7 ...
           + 2.6833087e1 * T.^8 ...
           - 4.5875709e0 * T.^9) * 10;
  end

  % Osborne, Abraham, Weinstock, 1951, T = 0.5 .. 1.5 K
  function P = pmelt_osb(T)
    P = (26.8e0 + 13.1e0 * T.^2) * 1.01325e0;
  end

  % Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955)) T = 2 - 31 K
  function P = pmelt_mills(T)
    P = (25.16e0 + 20.08201e0 * T.^1.517083e0) * 0.980665e0;
  end


  % plot the data:
  t1=9e-4:0.001:0.25;
  plot(t1, pmelt_gr(t1)+DGR, 'c-', 'linewidth', 3);

  t2=9e-4:0.001:1.0;
  plot(t2, pmelt_plts(t2), 'r-', 'linewidth', 3);

  t3=0.5:0.01:1.5;
  plot(t3, pmelt_osb(t3), 'g-', 'linewidth', 3);
  plot(t3, pmelt_osb(t3)+DGR, 'g-', 'linewidth', 3);

  t4=2:0.1:3.1;
  plot(t4, pmelt_mills(t4), 'm-', 'linewidth', 3);
  plot(t4, pmelt_mills(t4)+DGR, 'm-', 'linewidth', 3);

  % do PLTS - Mills interpolation
  t1 = 0.9:0.01:1; t2 = 1.4:0.01:1.5; t3=2:0.01:3;
  tint = [t1 t2 t3];
  pint = [pmelt_plts(t1) pmelt_osb(t2)  pmelt_mills(t3)];
  p1 = polyfit(tint, pint, 3)

  % plot it
%  t5=0.9:0.01:3;
%  plot(t5, polyval(p1, t5), 'b-', 'linewidth', 1);

  % define smooth functions for PLTS and greywall scales:

  function P = pmelt_ext_plts(T)
    T1 = 1.0;
    T2 = 1.1;
    T3 = 2.0;
    T4 = 3.0;
    pi = [-2.8399 21.1585 -3.5740 25.1894];
    P=zeros(size(T));
    ii = find(T<=T1);          P(ii) = pmelt_plts(T(ii));
    ii = find(T>=T2 & T<=T3 ); P(ii) = polyval(pi, T(ii));
    ii = find(T>=T4);          P(ii) = pmelt_mills(T(ii));
    ii = find(T>T1 & T<T2);
      P(ii) = ((T(ii)-T1).*polyval(pi, T(ii)) ...
             + (T2-T(ii)).*pmelt_plts(T(ii)))/(T2-T1);
    ii = find(T>T3 & T<T4);
      P(ii) = ((T(ii)-T3).*pmelt_mills(T(ii)) ...
             + (T4-T(ii)).*polyval(pi, T(ii)))/(T4-T3);
  end

  function P = pmelt_ext_gr(T)
    T1 = 0.25;
    T2 = 0.27;
    P=zeros(size(T));
    ii = find(T<=T1); P(ii) = pmelt_gr(T(ii));
    ii = find(T>=T2); P(ii) = pmelt_ext_plts(T(ii));
    ii = find(T>T1 & T<T2);
      P(ii) = ((T(ii)-T1).*pmelt_plts(T(ii)) ...
             + (T2-T(ii)).*pmelt_gr(T(ii)))/(T2-T1);
  end


  tt=0:0.001:3.1;
  plot(tt, pmelt_ext_plts(tt), 'k-');
  plot(tt, pmelt_ext_gr(tt)+DGR, 'k-');

  % test library functions
  plot(tt, he3_pmelt_plts(tt), 'b--');
  plot(tt, he3_pmelt(tt)+DGR, 'b--');

  xlabel('temperature, K');
  ylabel('p, bar');



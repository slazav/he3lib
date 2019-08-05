#!/usr/bin/octave -qf

  % interpolation of the melting pressure using various datasets

  addpath ../../matlab
  figure; clf; hold on;

  DGR = 0.1; % shift Greywall data

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

  % plot the data:
  t1=0.9:0.01:2;
  plot(t1, pmelt_gr(t1/1e3)+DGR, 'b-', 'linewidth', 3);
  plot(t1, pmelt_plts(t1/1e3), 'r-', 'linewidth', 3);

  t2=0:0.01:2;
  plot(t2, pmelt_gr(t2/1e3)+DGR, 'b-', 'linewidth', 1);
  plot(t2, pmelt_plts(t2/1e3), 'r-', 'linewidth', 1);

  tf = 0.9:0.1:2;
  pp1 = polyfit(tf, pmelt_gr(tf/1e3), 2)
  pp2 = polyfit(tf, pmelt_plts(tf/1e3), 2)

  plot(t2, polyval(pp1, t2)+DGR, 'k-');
  plot(t2, polyval(pp2, t2), 'k-');

  ylim([34.2 34.6])

  xlabel('temperature, K');
  ylabel('p, bar');


function plot_pars()
  find_figure('He3 F0'); clf;
  press=0:0.1:35;

  % data from VW (Wheately 1975)
  pdat = [0   3   6   9   12  15  18  21  24  27  30  33  34.39];
  fdat = [695 723 733 742 747 753 757 755 756 755 754 755 753] * -1e-3;

  %parameter F_0^a, following D.Einzel JLTP 84 
  pdatE=[0 10 20 30];
  fdatE=-[0.7 0.74 0.76 0.75];

  %
  hold on; title('F0a');
  plot(press, he3_f0a(press), 'b-');
  plot(press, interp1(pdatE,fdatE, press,'spline'), 'r-');
  plot(pdat,  fdat, 'g*');
  plot(press, F0a1(press), 'm-');

  % http://spindry.phys.northwestern.edu/3He Calculator/Tabfoot.html
  % Ramm, H., P. Pedroni, J.R. Thompson, and H. Meyer, J. Low Temp. Phys. 2, 539 (1970).
  % The fit coefficients are from Halperin, W.P., and E. Varoquax, Helium Three ed.
  %   W.P. Halperin, and L.P. Pitaevskii, Elsevier (1990). 
  p = [-1.823e-6 2.057e-4 -6.232e-3 -7.007e-1];
  plot(press, polyval(p, press), 'k-')

  plot(pdatE,  fdatE, 'r*');

  legend('Greywall-1983, Mukharskii program', ...
         'Einzel JLTP 84 -- Samuli', ...
         'Vollhardt-Wolfle - Wheately-1975', ...
         'Erkki, Thuneberg -- ve', ...
         'Ramm-1970, Halperin-1990');

end


function res = molvol(p)
% molar volume of 3He in cm^3 versus pressure, bar 
% from Erkki's calculator 
  res = 36.837231+p.*(-1.1803474+p.*(0.083421417+p.*(-0.0038859562+p.*(0.00009475978-p.*0.00000091253577))));
end

function res = F0a1(p)
  % fermi parameter F0a as a function of pressure, Thuneberg01 fit 
  res = -0.909+0.0055*molvol(p);
end

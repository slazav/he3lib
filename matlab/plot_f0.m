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

  plot(pdatE,  fdatE, 'r*');

  legend('Greywall 1983, Mukharskii program', ...
         'Einzel JLTP 84 -- Samuli', ...
         'Vollhardt-Wolfle - Wheately(1975)', ...
         'Erkki, Thuneberg -- ve')

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

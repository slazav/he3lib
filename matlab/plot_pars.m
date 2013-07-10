function plot_pars()
  find_figure('He3 parameters'); clf;
  press=0:0.1:35;

  subplot(3,2,1); hold on; title('molar volume')
  plot(press, he3_vm(press), 'b-');
  xlabel('pressure, bar');
  ylabel('molar volume, cm**3/mole');
  xlim([0 35]);

  subplot(3,2,2); hold on; title('atom mass and effective mass')
  plot(press, he3_meff(press), 'b-');
  plot([0 35], he3_amass*[1 1], 'r-');
  xlabel('pressure, bar');
  ylabel('mass, g');
  xlim([0 35]);

  subplot(3,2,3); hold on; title('Fermi momentum')
  plot(press, he3_pf(press), 'b-');
  xlabel('pressure, bar');
  ylabel('Fermi momentum, g*cm/s');
  xlim([0 35]);

  subplot(3,2,4); hold on; title('Fermi velocity')
  plot(press, he3_vf(press), 'b-');
  xlabel('pressure, bar');
  ylabel('Fermi velocity, cm/s');
  xlim([0 35]);

  subplot(3,2,5); hold on; title('R-Gas constant')
  plot(press, he3_gammaf(press), 'b-');
  xlabel('pressure, bar');
  ylabel('gamma, 1/(K*mol)');
  xlim([0 35]);

  subplot(3,2,6); hold on; title('Density of state')
  plot(press, he3_dnde(press), 'b-');
  xlabel('pressure, bar');
  ylabel('dnde, 1/sgs??');
  xlim([0 35]);

end
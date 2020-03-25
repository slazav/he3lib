#!/usr/bin/octave-cli

  ttc=0.1:0.1:1;
  p=10;

  nu_b = he3_nu_b(ttc,p);

  fprintf('ttc nu_b,kHz\n');
  fprintf('%0.1f %7.3f\n', [ttc; nu_b/1000]);

  fprintf('\n');
  fprintf('m_3 =%.3e\n', he3_amass);

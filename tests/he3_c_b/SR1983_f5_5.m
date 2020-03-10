function SR1983_f5_5
  addpath ~/PROG/he3lib/octave
  figure; clf; hold on;

  plot1(1.8, 17.25, 'b');
  plot1(1.6, 4.68,  'r');


  xlim([0 1])
  ylim([0 4])
  xlabel('T/Tc')
  ylabel('C(T)/C_N(Tc) / (T/Tc)^2')

  legend('\Delta C/C_N = 1.8', '', '1.6', '',
         'location', 'southeast')
  legend('boxoff')
  print -deps -color "-S400,350" SR1983_f5_5.eps

end


function plot1(fn, p, c)
  ttc = 0.05:0.01:1;
  f = fopen(['SR1983.f5.5_', sprintf('%.1f', fn), '.txt']);
  x = textscan(f, '%f %f');
  gap = he3_trivgap(ttc,p);
  dgap2 = he3_trivdgap2(ttc,p);
  plot(ttc, he3_yosida_c(ttc, gap, dgap2)./ttc.^2, [c '-']);
  plot(x{1}, x{2}, [c '*']);
end



function greywall86_f18
  addpath ~/PROG/he3lib/octave
  figure; clf; hold on;

  plot1(7, 00.00, 'r');
%  plot1(8, 02.18, 'g');
%  plot1(9, 05.21, 'b');
  plot1(4, 10.25, 'g');
%  plot1(5, 14.95, 'g');
  plot1(6, 20.30, 'b');
% plot1(1, 25.31, 'r');
% plot1(2, 29.08, 'g');
  plot1(3, 33.95, 'm');

  xlim([0.2,1])
  ylim([1,4])

  xlabel('T/Tc')
  ylabel('C(T)/C_N(Tc) / (T/Tc)^2')

  legend('0 bar', '', '10 bar', '', '20 bar', '', '33 bar', '',
         'location', 'southeast')
  legend('boxoff')
  print -deps -color "-S400,350" greywall86_f18.eps

end


function plot1(fn, p, c)
  ttc = 0.05:0.01:1;
  f = fopen(['greywall86_f18_', num2str(fn), '.txt']);
  x = textscan(f, '%f %f');
  gap = he3_trivgap(ttc,p);
  dgap2 = he3_trivdgap2(ttc,p);
  plot(ttc, he3_yosida_c(ttc, gap, dgap2)./ttc.^2, [c '-']);
  plot(0.2+0.8*x{1}, x{2}, [c '*']);
end

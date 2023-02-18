#!/usr/bin/octave-cli -f
pkg load he3lib

figure; hold on;

ttc=0.1:0.001:1;
p=0;
gap=he3_gap(ttc, p);

y0 = he3_yosida(ttc, gap, 0);
y1 = he3_yosida(ttc, gap, 1);
y2 = he3_yosida(ttc, gap, 2);
y3 = he3_yosida(ttc, gap, 3);

y5 = (8D0/15.0*y2./y1 + 5/8.0*y3./y2)./y2;
y6 = sqrt(y2.*y0);

ts = ttc.*(0.9074 - 0.0075*ttc - 0.0216*ttc.^2 + 0.1396*ttc.^3 - 0.0611*ttc.^4);

y0a = ttc
ii=ttc<0.94;
y0a(ii)  = exp(-1.76388./ts(ii)).*(3.454 - 0.88*ts(ii) + 4.625*ts(ii).^2 - 1.367*ts(ii).^3)./sqrt(ts(ii));
ii=(ttc>=0.94);
y0a(ii) = 1.985*ts(ii) - 0.985;

y5a = ttc
ii=ttc<0.80;
y5a(ii) = exp(1.76388./ts(ii)).*(0.10177 + 1.1958*ts(ii) - 1.425*ts(ii).^2 + 0.392*ts(ii).^3)./sqrt(ts(ii));
ii=(ttc>=0.80);
y5a(ii) = exp(1.76388./ts(ii)).*(0.19847 + 0.335*sqrt(1-ts(ii)))./sqrt(ts(ii))

y6a = ttc
ii=ttc<0.90;
y6a(ii) = exp(-1.76388./ts(ii)).*(2.402 + 0.4467*ts(ii) - 2.117*ts(ii).^2 + 4.1*ts(ii).^3)
ii=(ttc>=0.90);
y6a(ii) = 1.0 - sqrt(1.0-ts(ii)).*(-4.517 + 13.275*ts(ii) - 7.5*ts(ii).^2)

plot(ttc, y0a./y0, 'r-')
plot(ttc, y5a./y5, 'b-')
plot(ttc, y6a./y6, 'b-')

print -dpng he3b_tau.png "-S640,480"

%
%
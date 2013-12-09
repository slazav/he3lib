#!/usr/bin/octave -qf

addpath ../matlab
% values from paper:
[bs, as] = textread('2009_surovtsev.dat', '%f %f',...
                  'commentstyle', 'shell');
bs=bs*180/pi;
ab=mean(as(2:4)./bs(2:4))

%%%
  b=0:0.01:180;
  i1=find(cosd(b)>-1/4);
  i2=find(cosd(b)<=-1/4);
  b1=b(i1);
  b2=b(i2);

  a11(i1)=sqrt( (1-cosd(b1)) .* (1+4*cosd(b1)) )/5;
  a11(i2)=zeros(size(i2));

  a12(i1)= (1-cosd(b1)) * sqrt(10)/25;
  a12(i2)= abs(4*cosd(b2).^2+31*cosd(b2)+15)/15/sqrt(10);

  a13(i1)= a11(i1) * sqrt(10)/5;
  a13(i2)= -sind(b2).*(1+4*cosd(b2)) * sqrt(10)/225;

  a1213=sqrt(a12.^2 + a13.^2);

  a23(i1) = (1-cosd(b1))*3/20;
  a23(i2) = abs(4*cosd(b2).^2+31*cosd(b2)+15)/40;

  a=max([a11; a1213; a23]);

%%%

ab = pi/180 / sqrt(10)


figure; clf; hold on;
plot(b,a11, 'r-');
plot(b,a1213, 'g-');
plot(b,a23, 'b-');
plot(b,a, 'r-', 'linewidth', 2);
plot(bs,as, 'b*');

plot([0 40], ab*[0 40], 'b-');
xlabel('\beta, deg');
xlabel('a(\beta)');
legend('a11', 'a12-13', 'a23', 'a', 'picture from paper', 'location', 'southeast')

t=0.12:0.002:0.26;
p=20;
f=250000;

D = he3_diff_perp_xx(t,p,f);
o0=2*pi*f;
ob=2*pi*he3_nu_b(t,p);

cperp=he3_swvel_per(p,t*he3_tc(p));
cpar=he3_swvel_par(p,t*he3_tc(p));
mu=0.25; % 1-cperp.^2./cpar.^2

figure; clf; hold on;
semilogy(t, (D*o0^3)./(2*ob.^2.*cpar.^2.*mu*ab)  , 'b.-');
xlabel('T/T_c')
ylabel('\beta_{max}, deg')


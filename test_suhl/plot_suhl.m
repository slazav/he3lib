#!/usr/bin/octave -qf

addpath ../matlab
[b, a] = textread('2009_surovtsev.dat', '%f %f',...
                  'commentstyle', 'shell');
b=b*180/pi;

ab=mean(a(2:4)./b(2:4))

figure; clf; hold on;
plot(b,a, 'b.-');
plot([0 40], ab*[0 40], 'r-');
xlabel('\beta, deg');
xlabel('a(\beta)');

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
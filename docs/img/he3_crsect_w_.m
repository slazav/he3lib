#!/usr/bin/octave-cli -qf

pkg load he3lib

figure; clf; hold on;
xlabel('P');
ylabel('<W>');



P=0:30;
v1=he3_crsect_w(P);

pp=[0,10,20,30];
v2=[66.6, 105.9, 123.7, 126.5];
v3=[74.3, 109.6, 118.6, 121.2];


subplot(2,2,1); hold on
plot(P,v1, 'm-');
plot(pp,v2, 'b*-');
plot(pp,v3, 'r*-');

xlabel('P, bar')
legend('he3\_crsect\_w',...
       '<W> Einzel-1984, tab3, sp-wave',...
       '<W>Einzel-1984, tab3, P&W',...
 'location', 'southeast');

subplot(2,2,2); hold on
g0a=[0.119, 0.108, 0.097, 0.094];
g0b=[0.135, 0.156, 0.139, 0.089];
plot(P,he3_scatt_g0(P), 'm-');
plot(pp,g0a, 'b*-');
plot(pp,g0b, 'r*-');

xlabel('P, bar')
legend('he3\_scatt\_g0',...
       '<W> Einzel-1984, tab3, sp-wave',...
       '<W>Einzel-1984, tab3, P&W',...
 'location', 'southwest');

subplot(2,2,3); hold on
d0a=[0.295, 0.301, 0.303, 0.301];
d0b=[0.268, 0.296, 0.285, 0.252];
plot(P,he3_scatt_d0(P), 'm-');
plot(pp,d0a, 'b*-');
plot(pp,d0b, 'r*-');

xlabel('P, bar')
legend('he3\_scatt\_d0',...
       '<W> Einzel-1984, tab3, sp-wave',...
       '<W>Einzel-1984, tab3, P&W',...
 'location', 'southwest');

subplot(2,2,4); hold on
l1a=[1.25 0.966 0.910 0.928];
l1b=[1.19 1.33, 1.33, 1.26];
plot(P,he3_scatt_l1a(P), 'm-');
plot(pp,l1a, 'b*-');
plot(pp,l1b, 'r*-');

xlabel('P, bar')
legend('he3\_scatt\_l1a',...
       '<W> Einzel-1984, tab3, sp-wave',...
       '<W>Einzel-1984, tab3, P&W',...
 'location', 'east');


#subplot(2,2,4); hold on
#plot(P,he3_crsect_wi(P), '-');
#plot(P,he3_crsect_wd(P), '-');
#plot(P,he3_crsect_wl(P), '-');
#xlabel('P, bar')
#legend('he3\_crsect\_wi', 'he3\_crsect\_wd', 'he3\_crsect\_wl',...
# 'location', 'northwest');



print he3_crsect_w.png -dpng "-S800,600" "-F:6"

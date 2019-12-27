%Demag simulation
%he3 0bar
%

% heat leak to helium
qh=550; %pW
qh=qh*1.0e-12;
qhd=qh;

k=1.18;  %Korringa in sK

cu=17; %copper volume in cc
nc=cu/7.11; %moles
cn0=3.2e-6*nc; % heat cap = cn0 B^2/T^2

a=200;  %sinter area in m^2

he=35;  %helium volume in cc
nh=he/36.8; %moles at 0 bar
ch0=nh*24;   %heat cap = ch0 T

%demag  profile
b0=7; %initial field
b1=0.5; %mid f1eld
bf=0.056;% final field

tt=30; % total time in hours
td=12; %demag time, to mid field in hours
tdd= 12; % demag time, mid to final field
sd=td*3600;
sdd=tdd*3600;
ttt=tt*3600; %secs
step=.5;

s=[0:step:(sd+sdd)]; % demag array
t=[0:step:ttt];% total seconds array
b=ones(size(t))*bf;
bstep=zeros(size(t));
demend1=sd/step; %end mid demag
demend=(sdd+sd)/step;  %end demag

for i= 1:demend1-1
    bstep(i)=step*(b1-b0)/sd;
    b(i)=b0+bstep(i)*(i-1);
end           %field array start
for i= demend1-1: demend-1
    bstep(i)=step*(bf-b1)/sdd;
    b(i)=b1+bstep(i)*(-demend1+i-1);
end           %field array

%start temperature
t0=0.006;  %K

w=zeros*(size(t)); %array for demag energy
tn=t0*ones(size(t));
te=t0*ones(size(t));
th=t0*ones(size(t));  %initial setup
ch=ch0*ones(size(t));
cn=cn0*ones(size(t));
koroverkap=ones(size(t));
boverth=ones(size(t));
b(1)=b0;
kor=zeros(size(t));
kap=zeros(size(t));

for i=2:length(t)
   
   cn(i)=cn0*b(i)^2/tn(i)^2;
   kor(i)=tn(i-1)*(te(i-1)-tn(i-1))/k; %cn(i) left out
   tn(i)= tn(i-1)+kor(i)*step + bstep(i-1)*tn(i-1)/b(i-1); %new tn
    
   kap(i)=a*(0.15*tn(i-1)+4940*tn(i-1)^3)/741; %kap link
   te(i)=(kap(i)*th(i-1)+cn(i-1)*tn(i-1)^2/k)/(kap(i)+cn(i-1)*tn(i-1)/k);
   %new te
   ch(i)=ch0*th(i-1);
   th(i)=th(i-1)+(qhd-kap(i)*(th(i-1)-te(i-1)))*step/ch(i);
   %new th
   
    koroverkap(i)= (kor(i)*cn(i))/(kap(i)*(th(i)-te(i)));
   boverth(i)= b(i)/th(i);
end  
     
printfile5











%Demag simulation
%he3 0bar
%

% heat leak to helium
qh=300; %pW
qh=qh*1.0e-12;
qhd=qh;

qc=300; % steady pW to copper
qc=qc*1e-12;
qcb=12;  % wick for eddy current prop Bdot^2

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

tt=72; % total time in hours
td=12; %demag time, to mid field in hours
wait=0; %wait time at mid field
tdd=6; % demag time, mid to final field
sd=td*3600;
sdd=tdd*3600;
ttt=tt*3600; %secs
step=1;

s=[0:step:(sd+sdd)]; % demag array
t=[0:step:ttt];% total seconds array
b=ones(size(t))*bf;
bstep=zeros(size(t));
demend1=sd/step; %end mid demag
w1=wait*3600/step; %wait steps
demend=(sdd+sd+wait*3600)/step;  %end demag

for i= 1:demend1-1
    bstep(i)=step*(b1-b0)/sd;
    b(i)=b0+bstep(i)*(i-1);
end           %field array start
for i=demend1-1:demend1-1+w1
    bstep(i)=0;
    b(i)=b1;
end
for i=demend1-1+w1: demend-1
    bstep(i)=step*(bf-b1)/sdd;
    b(i)=b1+bstep(i)*(-demend1-w1+i-1);
end           %field array

t0=0.006;   %K   %start temperature

w=zeros*(size(t)); %array for demag energy
tn=t0*ones(size(t));
te=t0*ones(size(t));
th=t0*ones(size(t));  %initial setup
ch=ch0*ones(size(t));
cn=cn0*ones(size(t));
koroverkap=ones(size(t));
boverth=ones(size(t));
b(1)=b0;
kap=zeros(size(t));
qcc=zeros(size(t));

s0=-cn0*0.5*b0^2/t0^2 + ch0*t0    %initial entropy
s=zeros(size(t));

for i=2:length(t)
   
   kor(i)=tn(i-1)*(te(i-1)-tn(i-1))/k; %cn(i) left out
   tn(i)= tn(i-1)+(kor(i))*step + bstep(i-1)*tn(i-1)/b(i-1); %new tn
   cn(i)=cn0*b(i-1)^2/tn(i-1)^2;
   kap(i)=a*(0.15*0.5*(th(i-1)^2-te(i-1)^2)+4940*0.25*(th(i-1)^4-te(i-1)^4))/741; %kap link
   ch(i)=ch0*th(i-1);
   th(i)=th(i-1)+(qhd-kap(i))*step/ch(i);  %new th
   
   %num=(a/741)*(0.0725*te(i-1)^2+0.75*4940*te(i-1)^4+0.0725*th(i-1)^2+0.25*4940*th(i-1)^4)+cn(i-1)*tn(i-1)^2/k + qc;
   %denom = (a/741)*(0.15*te(i-1)+4940*te(i-1)^3)+cn(i-1)*tn(i-1)/k;
   qcc(i)=qcb*bstep(i)^2/step^2;
   qc1=qc+qcc(i);
   %te(i)=num/denom;
   %new te
    %te(i)=(kap(i)+qc+ch(i-1)*th(i-1)^2/k)/(ch(i-1)*th(i-1)/k);
        
    %num= kor(i)*ch(i-1) - qc - kap(i);
    %dem= ch(i-1)*tn(i-1)/k + (a/741)*(0.15*te(i-1)+4940*te(i-1)^4);
    %te(i)=te(i-1)-num/dem;
    num=(a/741)*(0.15*th(i-1)^2+4940*th(i-1)^4)+cn(i-1)*tn(i-1)^2/k + qc1;
    denom = (a/741)*(0.15*th(i-1)+4940*th(i-1)^3)+cn(i-1)*tn(i-1)/k;
   te(i)=num/denom;
   
   s(i)=(th(i)-tn(i))*step*kap(i)/tn(i)^2 + qhd*step/th(i)+qc1*step/te(i);
   
   koroverkap(i)= (kor(i)*cn(i))/(kap(i)+qc1);
   boverth(i)= b(i)/th(i);
end
slost=0;
for i=2:length(t)
    slost=slost+s(i);
end
slost
printfile5











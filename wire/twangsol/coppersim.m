%Platinum NMR simulation
%31.25 kHz (3.415 mT)
%

%scale heat for 31.25kHz
qs=4.233; %pW
qs=qs*1.0e-12;

%Korringa in sK
k=1.18;

%alpha
alpha=0.0041;

%steady heat leak
q0=4.0;  %pW

%times in hours
step=.002;
t=[0:step:6.5];
fudge=17; %thermometer slowing

%heat leak total
q=q0*ones(size(t));
for i=200:length(t)
   if t(i)>20.49
      if t(i)<20.95   %20.873
         q(i)=q(i)+23.6;
      end
   end
   if t(i)>24.556
      if t(i)<25.14
         q(i)=q(i)+.73;
      end
   end
     
  if t(i)>25.206
     if t(i)<27.69
         q(i) = q(i)+3.83;
      end
   end
end

q=q*1.0e-12;  %scale to watts

%initialise
step=step*3600;
t=t*3600;   %seconds
chi=zeros(size(t));
invtn=zeros(size(t));
invte=zeros(size(t));

%start values
tn=47.0e-6;    %kelvin

invtn(1)=1/tn;
chi(1)=20;
invte(1)=invtn(1)*qs/(qs+q(1));

for i=2:length(t)
   invtn(i)=invtn(i-1)-step*q(i-1)/(k*qs);
   invte(i)=invtn(i)*qs/(qs+q(i-1));
   chi(i)=chi(i-1)-(step/fudge)*(-1+chi(i-1)/(invte(i)*alpha));
end

chi=chi+8;  %background
invte=alpha*invte;
invtn=alpha*invtn;
printfile4
%Platinum NMR simulation
%31.25 kHz (3.415 mT)
%


%scale heat
qs=4.233; %pW
%qs=qs*10^(-12)

%Korringa
k=1.18;

%alpha
alpha=0.0043;

%steady heat leak
q0=2.0;  %pW

%times in hours
t=[0:.01:50];

%initialise
t=t*3600;   %seconds
chi=zeros(length(t));
invtn=zeros(length(t));
invte=zeros(length(t));

%start values
tn=12e-6;    %kelvin

invtn(1)=1/tn;
chi(1)=40;
invte(1)=chi(1);

chi(1)
chi(2)

end

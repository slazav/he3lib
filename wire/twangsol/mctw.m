% dilute saturated solution calibration
%
% enter wire-specific data here and type 
%
% mccal(f)
%
% in command window, where f is desired frequency
% 
% To modify temperature points, edit file temps.m
%
function mccal(f)
global  f a alpha rhorat rho

%enter data
%frequency
%f=600; entered from command window

%tantalum wire
rhow=16.7;
diam=124;  %microns
alpha=2.2;

a=diam/2e6; %radius in m

%helium data for saturated solution
conc=6.65;
vol=27.58;
mstar=2.46;
vol3=vol*(1+.0028*conc);
rho=conc*3*mstar*.01/vol3;

rhorat=rho/rhow;

%get temperatures
T=[2,8,15,35,65,100]

T=T/1000; %convert to Kelvin
n=length(T);

%pick up Stokes G and mfp L as arrays

for c=1:n,
   t=T(c);
   G(c)=a/pend(t);
   L(c)=mfp(t);
	visc(t)*t^2
end

%work out shift and width

for c=1:n,
   g=G(c);
   l=L(c);
   F2(c)= f*rhorat*stokesk3(g,l)*(1-1.14*rhorat*stokesk2(g,l));
   F1(c)= f*rhorat*0.5*stokesk2(g,l)*(1-0.75*rhorat*stokesk2(g,l));
   Ft(c)= f*rhorat*stokesk1(g);
   %Ftt(c)=f*rhorat*stokesk3(g,l);
   
end

F2

Ft




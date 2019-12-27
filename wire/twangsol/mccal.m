% dilute saturated solution calibration
%
% enter wire-specific data here and type 
% mccal
% in command window, where f is desired frequency
% 
% To modify temperature points, edit file temps2.m
%


global  f a alpha rhorat rho

%enter data
%frequency
f=656; %not entered from command window

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
temps2

T=T/1000; %convert to Kelvin
n=length(T);

%pick up Stokes G and mfp L as arrays

for c=1:n,
   t=T(c);
   G(c)=a/pend(t);
   L(c)=mfp(t);
end

%work out shift and width
f
for c=1:n,
   g=G(c);
   l=L(c);
   F2(c)= f*rhorat*stokesk3(g,l)*(1-1.14*rhorat*stokesk2(g,l));
   F1(c)= f*rhorat*0.5*stokesk2(g,l)*(1-0.75*rhorat*stokesk2(g,l));
   %Ft(c)= f*rhorat*stokesk1(g);
   %Ftt(c)=f*rhorat*stokesk3(g,l);
   Inv(c)=1/T(c);
end

P=polyfit(F2,Inv,5); %polynomial fit to data 
PP=fliplr(P);

poly_coeffs_InverseT_vs_width =PP

for c=5:250,
   df=c/10;
   DF=[1,df,df^2,df^3,df^4,df^5];
   Inverset(c)=DF*PP';
   Width(c)=df;
end

plot (Inv,F2,'r+',Inverset,Width, 'b-')



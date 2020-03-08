% dilute saturated solution calibration
clear
global  f a alpha rhorat rho ee

addpath twangsol
addpath ../octave

% input data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = 'ian MC TW';
rhow=16.7; % tantalum density in g/cc 16.7, quartz 2.659, NbTi 6.05
diam=140; % diameter in microns 124,  fork 101,  micro 13.5, triple 4.5
f=1092;  %frequency in Hz  385

%Temp=[5:.5:10,10.5:1:30];
T=[2:.1:5,5.2:.2:10,10.5:.5:19.5, 20:1:90];% temperatures required in mK
T=fliplr(T);
file = 'new/MC_calibration.dat';  %file to store data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=2.2;   % mfp fudge.
ee=10;  % switch strength to ballistic
a=diam/2e6;  %radius in m
%helium data for saturated solution
%conc=6.65;
vol=27.58;
mstar=2.46;
T=T/1000; %convert to Kelvin
n=length(T);
%pick up Stokes G and mfp L as arrays
for c=1:n,
   t=T(c);
   C(c)=0.066+0.5056*t^2-0.2488*t^3+18.22*t^4-74.22*t^5;
   vol3(c)=vol*(1+.286*C(c));
   rh(c)=C(c)*3.016*mstar/vol3(c);
   rhora(c)=rh(c)/rhow;
   rho=rh(c);
   G(c)=a/pend(t);
   cratio=0.0665/C(c);
   L(c)=visc(t)/105.5302*((cratio)^(4/3));
    
end
%work out shift and width
for c=1:n,
   rho=rh(c);
   rhorat=rhora(c);
   g=G(c);
   l=L(c);
   F2(c)= f*rhorat*stokesk3(g,l)*(1-1.14*rhorat*stokesk2(g,l));
   F1(c)= f*rhorat*0.5*stokesk2(g,l)*(1-0.75*rhorat*stokesk2(g,l));
   Inv(c)=1/T(c);
end
T=T*1000;  % back to mK

P=polyfit(F2,Inv,5); %polynomial fit to data 
PP=fliplr(P);
%poly_coeffs_InverseT_vs_width =PP
wmax=max(F2);
for c=1:wmax,
   df=c;
   DF=[1,df,df^2,df^3,df^4,df^5];
   Inverset(c)=DF*PP';
   Width(c)=df;
end

plot (Inv,F2,'r+',Inverset,Width, 'b-')

% save data
fid=fopen(file,'w')
fprintf(fid, '# %s\n', name);
fprintf(fid, '#  frequency %9.2f Hz  \n', f);
fprintf(fid, '#  density  %9.3f g/cc  \n', rhow);
fprintf(fid, '#  diameter %8.1f microns  \n\n', diam);

OUT=[T;F2;F1];
fprintf(fid,'#Temperature(mK)  Width        Shift\n');
fprintf(fid,'%9.4f     %9.4f     %9.4f \n', OUT);

fprintf(fid,'\n\n# 5th order polynomial to give 1/T(mK) as powers of width\n');
fprintf(fid,'%14.6e \n', PP);

fclose(fid)


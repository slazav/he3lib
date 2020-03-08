% normal 3He calibration , variable freq and CHH viscosity 
%(valid up to 50 mK?)

clear
global  f a alpha rhorat rho ee

addpath twangsol
addpath ../octave

% input data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = 'cylinder programme using wide line treatment, CHH and no slip fudge';
rhow=16.7;   % tantalum density in g/cc 16.7, quartz 2.659, NbTi 6.05
diam=124;     % diameter in microns 124,  fork 101,  micro 13.5, triple 4.5
fv=813;     % vacuum frequency in Hz 
p=0;          % pressure in bar

%Temp=[5:.5:10,10.5:1:30];
T=[.93, .95, .975,1:.1:5,5.2:.2:10, 11:1:19, 20:2:50];% temperatures required in mK
T=fliplr(T);
file = 'new/NormalHe3_twanger_low_exact.dat';  %file to store data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=1.9;
ee=1; %switch to ballistic
a=diam/2e6; %radius in m
%helium data
vol=he3_vm(p);
rho=3.016/vol;
rhorat=rho/rhow;
eovl=0.2*(6.023e29/vol)^(4/3)*(3*9.8696)^(1/3)*1.0546e-34;
T=T/1000; %convert to Kelvin
n=length(T);
tc=he3_tc(p)/1000;
%check always in normal state (optional)
for c=n-10:n
   T(c)=max([T(c) tc]);
end
%pick up Stokes G and mfp L as arrays
for c=1:n,
   t=T(c);
   vis=1/(visca(p)*t^2+viscb(p));  %CHH viscosity
    %shift using approx formula for fork; remove this section for wires  
 %  x=log10(t*1000);
  % F1a(c)=(1990-2165*x +1150*x^2 -247*x^3 +8.65*x^4+2.57*x^5);    
   %f=fv-F1a(c);
   f=fv;      
   pen(c)=sqrt(vis/(2000*pi*rho*f));
   g=a/pen(c);
   l=vis/eovl;
   f=fv*(sqrt(1+(rhorat*stokesk(g)/2)^2)-rhorat*stokesk(g)/2);
   %f=fv*(sqrt(1+(rhorat*stokesk2(g,l)/2)^2)-rhorat*stokesk2(g,l)/2);
   F1a(c)=f;
   
  for k=1:5
    pen(c)=sqrt(vis/(2000*pi*rho*f));
    g=a/pen(c);
    l=vis/eovl;
    f=fv*(sqrt(1+(rhorat*stokesk(g)/2)^2)-rhorat*stokesk(g)/2);
    %f=fv*(sqrt(1+(rhorat*stokesk2(g,l)/2)^2)-rhorat*stokesk2(g,l)/2);
  end
   F1(c)=f;  %correct frequency
   
   F2a(c)=rhorat*stokesk1(g)*f;
   sq1=f*sqrt((rhorat*(stokesk1(g)+stokesk(g)))^2 + 4);
   sq2=f*sqrt((rhorat*(stokesk1(g)-stokesk(g)))^2 + 4);
   
  %F2a(c)=rhorat*stokesk3(g,l)*f;
   %sq1=f*sqrt((rhorat*(stokesk3(g,l)+stokesk2(g,l)))^2 + 4);
   %sq2=f*sqrt((rhorat*(stokesk3(g,l)-stokesk2(g,l)))^2 + 4);
   
   F2(c)= F2a(c)+sq1 - sq2;  %exact width
   %F1(c)= f*rhorat*0.5*stokesk2(g,l)*(1-0.75*rhorat*stokesk2(g,l));
   Inv(c)=1/T(c);
end

T=T*1000;  % back to mK
P=polyfit(F2,Inv,7); %polynomial fit to data
PP=fliplr(P);
%poly_coeffs_InverseT_vs_width =PP
wmax=max(F2);
for c=1:wmax,
   df=c;
   DF=[1,df,df^2,df^3,df^4,df^5,df^6,df^7];
   Inverset(c)=DF*PP';
   Width(c)=df;
end

plot (Inv,F2,'r+',Inverset,Width, 'b-')
% save data
fid=fopen(file,'w')
fprintf(fid, '# %s\n', name);
fprintf(fid, '# vacuum frequency %9.2f Hz  \n', fv);
fprintf(fid, '#  density  %9.3f g/cc  \n', rhow);
fprintf(fid, '#  diameter %8.1f microns  \n\n', diam);

OUT=[T;F2;F1;F2a;F1a];
fprintf(fid,'#Temperature(mK)  Width    Frequency  approx width  approx freq\n');
fprintf(fid,'%9.3f    %9.3f     %9.3f    %9.3f   %9.3f \n', OUT);

fprintf(fid,'\n\n# 7th order polynomial to give 1/T(mK) as powers of width\n');
fprintf(fid,'#%14.6e \n', PP);

fclose(fid)






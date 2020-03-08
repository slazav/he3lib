% B-phase 3He calibration 
clear
global  f a alpha rhorat rho ee 

addpath twangsol
addpath ../octave

% input data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = 'triplemicro';
rhow=6.05; % tantalum density in g/cc 16.7, quartz 2.659, NbTi 6.05
diam=4.5; % diameter in microns 124,  fork 101,  micro 13.5, triple 4.5
f=260;  %frequency in Hz
p=0;   % pressure in bar

T=[1:-.002:0.14];% temperatures required in mK

file = 'new/SuperfluidHe3_cal.dat';  %file to store data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
al=1.9;  % value of alpha (modify if required)
a=diam/2e6; %radius in m
ee=1; % ballistic switch wick

tc    = he3_tc(p)/1000;   %in kelvin
vol   = he3_vm(p);
rho3  = he3_rho(p);
ff    = he3_f1s(p);  %F1s parameter

eovl=0.2*(6.023e29/vol)^(4/3)*(3*9.8696)^(1/3)*1.0546e-34;
vistc=1/(visca(p)*tc^2+viscb(p));

%get relative temperatures 
T=T/1000;   %kelvin
n=length(T);
for c=1:n       %ensure all below Tc
    if T(c)>tc; 
        T(c)=tc;
    end
end
TR=T/tc;    %reduced temperature
%pick up Stokes G and mfp L as arrays
for c=1:n,
    vis=vistc*redvis(TR(c));     %fudged effective viscosity from CHH
   y0(c)=yosidy0(TR(c));
   y5=yosidy5(TR(c));
   y6=yosidy6(TR(c));
   rho(c)=rho3*(1+ff/3)*y0(c)/(1+ff*y0(c)/3);    %effective density
   rhorat(c)=rho(c)/rhow;                        %density ratio
   pen(c)=sqrt(vis/(2000*pi*rho(c)*f));          %penetration depth
   zeta(c)=0.5*y5*vis/eovl;                      %  effective slip length
    L(c)=vis/(eovl*y6);                          %mean free path
   alp(c)=1.156*al/(y5*y6);                      %effective alpha
   G(c)=a/pen(c);                                %gamma for Stokes
%work out shift and width
   alpha=alp(c);
   z=zeta(c);
   g=G(c);
   l=L(c);
   F2(c)= f*rhorat(c)*stokesk3s(g,l,z);%*(1-1.14*rhorat(c)*stokesk2s(g,l,z));
   F1(c)= f*rhorat(c)*0.5*stokesk2s(g,l,z);%*(1-0.75*rhorat(c)*stokesk2s(g,l,z));
   Inv(c)=1/T(c)/1000;      %1/T in mK
end
LL=log10(F2);
plot(Inv,LL)

P=polyfit(LL,Inv,7); %polynomial fit to data 
PP=fliplr(P);
%poly_coeffs_InverseT_vs_logwidth =PP
w=max(LL)-min(LL);
wlow=min(LL);
for c=1:200,
    df=wlow +c*w/200;
    DF=[1,df,df^2,df^3,df^4,df^5,df^6,df^7];
   Inverset(c)=DF*PP';
   LogW(c)=df;
end
plot (Inv,LL,'r+',Inverset,LogW, 'b-')
T=T*1000;  % back to mK
% save data
fid=fopen(file,'w')
fprintf(fid, '# %s\n', name);
fprintf(fid, '#  frequency %9.2f Hz  \n', f);
fprintf(fid, '#  density  %9.3f g/cc  \n', rhow);
fprintf(fid, '#  diameter %8.1f microns  \n\n', diam);

OUT=[T;F2;F1];
fprintf(fid,'# Temperature(mK)  Width        Shift\n');
fprintf(fid,'%9.4f     %9.4f     %9.4f \n', OUT);

fprintf(fid,'\n\n # 7th order polynomial to give 1/T(mK) as powers of log10width\n');
fprintf(fid,'# %14.6e \n', PP);

fclose(fid)

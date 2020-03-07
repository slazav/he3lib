% normal 3He calibration , variable freq and deepak or CHH viscosity 
%(CHH valid up to 50 mK?)  Cylinder approx with slip

clear
global  f a alpha rhorat rho ee

% input data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = 'TF 0 bar using wide line treatment and CHH viscosity';
rhow=2.659;   % tantalum density in g/cc 16.7, quartz 2.659, NbTi 6.05
diam=101;     % diameter in microns 124,  fork 101,  micro 13.5, triple 4.5
fv=32700;     % vacuum frequency in Hz 
p=0;          % pressure in bar

%Temp=[5:.5:10,10.5:1:30];
T=[.93, .95, .975,1:.1:5,5.2:.2:10, 10.2:.2:19, 20:1:30]%100,110:10:1600];% temperatures required in mK
T=fliplr(T);
file = 'new/NormalHe3_twanger_high_exact.dat';  %file to store data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=1.9;   %1.9
ee=10; %switch to ballistic   %1
a=diam/2e6; %radius in m
%helium data
vol=volume(p);
rho=3.016/vol;
rhorat=rho/rhow;
eovl=0.2*(6.023e29/vol)^(4/3)*(3*9.8696)^(1/3)*1.0546e-34;
T=T/1000; %convert to Kelvin
n=length(T);
tc=tche3(p)/1000;
%check always in normal state (optional)
for c=n-10:n
   T(c)=max([T(c) tc]);
end
%pick up Stokes G and mfp L as arrays
for c=1:n,
   t=T(c);
   vis=1/(visca(p)*t^2+viscb(p));  %CHH viscosity
   
   %vis=(1/t^2)*1e-7*(2.02717+4.7353*t+28.22503*t^2-24.7364*t^3+22.39929*t^4-6.92457*t^5);
   %Deepak's polynomial for viscosity
   
   f=fv;      
   pen(c)=sqrt(vis/(2000*pi*rho*f));
   g=a/pen(c);
   l=vis/eovl;
   f=fv*(sqrt(1+(rhorat*stokesk2(g,l)/2)^2)-rhorat*stokesk2(g,l)/2);
   F1a(c)=f;
   
  for k=1:10
    pen(c)=sqrt(vis/(2000*pi*rho*f));
    g=a/pen(c);
    l=vis/eovl;
    f=fv*(sqrt(1+(rhorat*stokesk2(g,l)/2)^2)-rhorat*stokesk2(g,l)/2);
  end
   F1(c)=f;  %correct frequency
   
   F2a(c)=rhorat*stokesk3(g,l)*f;
   sq1=f*sqrt((rhorat*(stokesk3(g,l)+stokesk2(g,l)))^2 + 4);
   sq2=f*sqrt((rhorat*(stokesk3(g,l)-stokesk2(g,l)))^2 + 4);
   
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
fprintf(fid,'# Temperature(mK)  Width    Frequency  approx width  approx freq\n');
fprintf(fid,'%9.3f    %9.3f     %9.3f    %9.3f   %9.3f \n', OUT);

fprintf(fid,'\n\n# 7th order polynomial to give 1/T(mK) as powers of width\n');
fprintf(fid,'%14.6e \n', PP);

fclose(fid)






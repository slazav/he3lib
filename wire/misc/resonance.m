% nonlinear resonance lines
%
% allow for velocity dependence
%
%output in velocity in mm/s

%enter data
f0=420.1;  %centre frequency

delta=.2;  %width for low velocities
k1v= 3000; %-12460  %non-linear param for shift - velocity
k2v=5;  %non-linear param for width - velocity
k1d=250; %non-li param shift - drive
k2d=.5;  % nonlin param for width - drive

wick=[.05 .1 .2 .5 1 2 5 10]; % drive parameters

f=[410:0.005:430];  % f sweep

n=length(f);
m=length(wick);
v(1)= 0;
v(2)=0;
for j=1:m,
   
   for i=3:n,
      vv=(v(i-1)+v(i-2))/2;
      d=delta + k2d*wick(j)+ k2v*vv;
      shift=k1d*wick(j)+k1v*vv^2;
      v(i)=wick(j)*d*f(i)^2/((f(i)^2-f0^2+shift)^2+d^2*f(i)^2);
   end

plot (f,v)


OUT=[f; v];

fid=fopen('d:\aerocell\calcs\test.dat','w'); %open file for write 
fclose(fid); %file flushed!
fid=fopen('d:\aerocell\calcs\test.dat','a') %open file for append 

if j-1,     
fprintf(fid,'frequency   velocity\n\n');
end
fprintf(fid,'%7.3f\n',wick(j));
fprintf(fid,'%9.5f %9.5f\n',OUT);
end
fclose(fid)

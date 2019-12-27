% dilute saturated solution calibration
%
% check
%


%get temperatures
T=[2:.2:100];

T=T/1000; %convert to Kelvin
n=length(T);


for c=1:n,
   t=T(c);
   V(c)=visc(t)*t^2*1e7;
end
X=[10:10:100];
Y=[.32, .345, .38, .415, .45, .485, .52, .58, .62, .685];


% X and Y are rough data taken from Zeegers et al
% graph of eta - t^2 vs t

plot (T*1000,V,'g-',X,Y,'r+')



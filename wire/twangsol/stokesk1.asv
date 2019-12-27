% stokes width function

function k1=stokesk1(g)

if g<.1
   x0=.02;
   x1=.1;
   y0=2390;
   y1=149.9;
elseif g<.2
   x0=.1;
   x1=.2;
   y0=149.9;
   y1=48.66;
elseif g<.5
   x0=.2;
   x1=.5;
   y0=48.66;
   y1=12.06;
elseif g<1
   x0=.5;
   x1=1;
   y0=12.06;
   y1=4.565;
elseif g<2
   x0=1;
   x1=2;
   y0=4.565;
   y1=1.876;
elseif g<5
   x0=2;
   x1=5;
   y0=1.876;
   y1=.641;
else k1=sqrt(8)/g + 2/g^2;
   return
end

x=log(g);
x0=log(x0);
x1=log(x1);
y0=log(y0);
y1=log(y1);

y=y0 + (x-x0)*(y1-y0)/(x1-x0);

k1=exp(y);

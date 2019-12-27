% stokes shift function

function k=stokesk(g)

if g<.2
   x0=.1;
   x1=.2;
   y0=46.71;
   y1=19.72;
elseif g<.5
   x0=.2;
   x1=.5;
   y0=19.72;
   y1=7.347;
elseif g<1
   x0=.5;
   x1=1;
   y0=7.347;
   y1=3.961;
elseif g<2
   x0=1;
   x1=2;
   y0=3.961;
   y1=2.438;
elseif g<5
   x0=2;
   x1=5;
   y0=2.438;
   y1=1.615;
else k=1+sqrt(8)/g;
  return
end

x=log(g);
x0=log(x0);
x1=log(x1);
y0=log(y0);
y1=log(y1);

y=y0 + (x-x0)*(y1-y0)/(x1-x0);

k=exp(y);


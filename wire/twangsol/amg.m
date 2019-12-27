function k2 = amg(g,l)
 
a=.101e-3/2;
k=stokes2k(g);
k1=stokes2k1(g);
alpha=1.9;
ee=1;

b=0.25*0.579*l/a;
b=b*(1+ee*alpha*l/a)/(1+ee*l/a);
num= k-1;
denom= (1+g^2*b*k1)^2+g^4*b^2*(k-1)^2;
k2=num/denom+1;
end


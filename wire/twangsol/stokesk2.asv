% shift function with mfp fudge

function k2=stokesk2(g,l)
global f a alpha rhorat ee
 
k1=stokes2k1(g);
k=stokes2k(g);
b=0.25*0.579*l/a;
b=b*(1+ee*alpha*l/a)/(1+ee*l/a);
num= k-1;
denom= (1+g^2*b*k1)^2+g^4*b^2*(k-1)^2;
k2 = 1 + num/denom;


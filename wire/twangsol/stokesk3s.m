% width function with mfp fudge

function k3=stokesk3s(g,l,z)
global f a alpha rhorat ee 

k1=stokes2k1(g);
k=stokes2k(g);
b=0.25*z/a;
b=b*(1+ee*alpha*l/a)/(1+ee*l/a);
num= k1+ g^2*b*((k-1)^2+k1^2);
denom= (1+g^2*b*k1)^2+g^4*b^2*(k-1)^2;
k3=num/denom;


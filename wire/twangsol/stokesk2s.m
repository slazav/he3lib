% shift function with mfp fudge

function k2=stokesk2s(g,l,z)
global f a alpha rhorat ee 

k1=math_stokes_kp(g);
k=math_stokes_k(g);
b=0.25*z/a;
b=b*(1+ee*alpha*l/a)/(1+ee*l/a);
num= k-1;
denom= (1+g^2*b*k1)^2+g^4*b^2*(k-1)^2;
k2 = 1 + num/denom;


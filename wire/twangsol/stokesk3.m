% width function with mfp fudge

function k3=stokesk3(g,l)
global f a alpha rhorat ee

k1=math_stokes_kp(g);
k=math_stokes_k(g);
b=0.25*0.579*l/a;
b=b*(1+ee*alpha*l/a)/(1+ee*l/a);
num= k1+ g^2*b*((k-1)^2+k1^2);
denom= (1+g^2*b*k1)^2+g^4*b^2*(k-1)^2;
k3=num/denom;


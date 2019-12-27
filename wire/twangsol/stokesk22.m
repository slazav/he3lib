% shift function with mfp fudge

function k2=stokesk22(g,l)
global f a alpha rhorat

b=0.25*0.579*l/a;
b=b*(1+10*alpha*l/a)/(1+10*l/a);
num= stokesk(g)-1;
denom= (1+g^2*b*stokesk1(g))^2....
   +g^4*b^2*(stokesk(g)-1)^2;

k2 = 1/(1+3*l/a) + num/denom;
%k2 = 1 + num/denom;

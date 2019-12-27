%viscous penetration depth

function d=pendsol(t);
global rho f 

d=sqrt(viscsol(t)/(2000*pi*rho*f));
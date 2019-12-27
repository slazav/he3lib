%viscous penetration depth

function d=pend(t);
d=sqrt(visc(t)/(2000*pi*rho*f));
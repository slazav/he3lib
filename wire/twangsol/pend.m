%viscous penetration depth

function d=pend(t)
global rho f

d=sqrt(visc(t)/(2000*pi*rho*f));
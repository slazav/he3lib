% mean free path for saturated solution

function l=mfpsol(t) 
global t conc
l=visc(t)/105.5302*((0.0665/conc)^(4/3));



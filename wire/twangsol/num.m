% bulk number function for given t and mu

function f=num(x);

global t  mu

f=x.^(1/2).*(exp((x-mu)/t)+1).^(-1);



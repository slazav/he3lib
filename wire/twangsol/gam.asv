%heat capacity gamma helium3

%   = C/RT

function gam=gam(p)

G(1)=0.27840464e1;
G(2)=0.69575243e-1;
G(3)=-0.14738303e-2;
G(4)=0.46153498e-4;
G(5)=-0.53785385e-6;
G(6)=0.0;
gam=0;
pr=1;
for c=1:6
    gam=gam+G(c)*pr;
    pr=pr*p;
end

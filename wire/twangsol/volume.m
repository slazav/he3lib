%molar volume of helium3

%    in cm^3 per mole

function vol=volume(p)

A(1)=36.837231;
A(2)=-0.11803474e1;
A(3)=0.83421417e-1;
A(4)=-0.38859562e-2;
A(5)=0.94759780e-4;
A(6)=-0.91253577e-6;
vol=0;
pr=1;
for c=1:6
    vol=vol+A(c)*pr;
    pr=pr*p;
end

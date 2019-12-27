%transition temperature of helium3

%    in mK

function tc=tche3(p)

 T(1)=0.92938375;
 T(2)=0.13867188;
 T(3)=-0.69302185e-2;
 T(4)=0.25685169e-3;
 T(5)=-0.57248644e-5;
 T(6)=0.53010918e-7;
tc=0;
pr=1;
for c=1:6
    tc=tc+T(c)*pr;
    pr=pr*p;
end

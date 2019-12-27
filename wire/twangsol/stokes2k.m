function k=stokes2k(g)
% Stokes k using fortran method

if g>3.0
         ak=2.828427/g;
         bk=.3535534/g^3;
         ck=.552427/g^5;
         dk=2.9637719/g^7;
         ek=20.4632383/g^9;
         k = 1+ak+bk-0.5/g^4+ck-dk+12.8875857/g^8-ek;
        
else
         gg=g/2;
         M0=gg^2-gg^6/12+gg^10/2880;
         M1=gg^2-gg^6/36+gg^10/14400;
         E0=gg^4/2-gg^8/144+gg^12/86400;
         E1=gg^4/4-gg^8/576+gg^12/518400;
         AL=.5772158 +log(gg);
         A=-(AL*M0)+pi/4*E0-0.5*M1+(gg^2-gg^6/6.545454+gg^10/1261.313);
         B=pi/4*M0+AL*E0-0.5*(1-E1)-(gg^4*0.75-gg^8/69.12+gg^12/35265.31);
         C=-(pi/4*M1)+AL*(1-E1)+(gg^4/2.666666-gg^8/276.48+gg^12/211591.84);
         D=-(AL*M1)-pi/4*(1-E1)+(gg^2-gg^6/19.636363+gg^10/6306.57);
         k=1 +2*(A*C + B*D)/(gg^2*(C*C + D*D));
         
end


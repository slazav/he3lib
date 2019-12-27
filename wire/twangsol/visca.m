%normal state viscosity

%     eta= 1/( AT^2 + B)\
%  eta in Pa s; T in mK originally, now K
%  using the data given in Carless, Hall and Hook; JLTP 50,583,83\
%  table 1 , page 593 as smoothed by AMG.\
%  The data, converted to a Greywall T scale\
%  by multiplication of T by 1.12.\

function aa=visca(p)

if p<1.28
    a=0.38-.007*(1.28-p)/1.18;
elseif p<4.65
    a=0.424-.044*(4.65-p)/3.37;
elseif p<9.89
    a=0.495-.071*(9.89-p)/5.24;
elseif p<19.89
    a=0.603-.108*(19.89-p)/10;
else a=0.603+0.107*(p-19.89)/9.45;
end  
aa=a*10*1.12*1.12;
aa=aa*1e6;  %convert back to T in K



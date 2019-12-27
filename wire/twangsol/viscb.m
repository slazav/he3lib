%normal state viscosity

%     eta= 1/( AT^2 + B)\
%
%  using the data given in Carless, Hall and Hook; JLTP 50,583,83\
%  table 1 , page 593 as smoothed by AMG.\
%  The data, converted to a Greywall T scale\
%  by multiplication of T by 1.12.\

function b=viscb(p)

if p<1.28
    b=0.06-0.02*(1.28-p)/1.18;
elseif p<4.65
     b=0.19-0.13*(4.65-p)/3.37;
elseif p<9.89
    b=0.43-0.24*(9.89-p)/5.24;
elseif p<19.89
    b=0.94-0.56*(19.89-p)/10;
else b=0.94+0.56*(p-19.89)/9.45;
end    
b=b*10;


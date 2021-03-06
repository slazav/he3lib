% yosida functions Y6
% CHH expression: t is reduced temperature (T/Tc)
% ts is their scaled temperature

function y6=yosidy6(t)

ts=t*(0.9074-0.0075*t - 0.0216*t*t +0.1396*t^3 - 0.0611*t^4);

if t<0.90
    y6=exp(-1.76388/ts)*(2.402+0.4467*ts-2.117*ts*ts+4.1*ts^3);
else
    y6=1- sqrt(1-ts)*(-4.517+13.275*ts-7.5*ts*ts);
    return

end

    
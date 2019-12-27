% yosida functions Y0
% CHH expression: t is reduced temperature (T/Tc)
% ts is their scaled temperature

function y0=yosidy0(t)

ts=t*(0.9074-0.0075*t - 0.0216*t*t +0.1396*t^3 - 0.0611*t^4);

if t<0.94
    y0=exp(-1.76388/ts)*(3.454-0.88*ts+4.625*ts*ts-1.367*ts^3)/sqrt(ts);
else
    y0=1.985*ts -0.985;
    return

end

    
% yosida functions Y5
% CHH expression: t is reduced temperature (T/Tc)
% ts is their scaled temperature

function y5=yosidy5(t)

ts=t*(0.9074-0.0075*t - 0.0216*t*t +0.1396*t^3 - 0.0611*t^4);

if t<0.80
    y5=exp(1.76388/ts)*(0.10177+1.1958*ts-1.425*ts*ts+0.392*ts^3)/sqrt(ts);
else
    y5=exp(1.76388/ts)*(0.19847+0.335*sqrt(1-ts))/sqrt(ts);
    return

end

    
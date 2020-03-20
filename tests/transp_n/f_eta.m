% function f_eta(lambda), VW eq.2.68
function f = f_eta(l)
  f = pi^2/12;

  for n=1:100
    df = l/4 .* (4*n+3) ./ ((2*n+1)*(n+1))^2 ./ ((2*n+1)*(n+1) - l);
    f = f + df;
    if max(df) < 1e-8; break; end
  end
end

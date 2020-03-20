% function f_k(lambda), VW eq.2.70
function f = f_k(l)
  f = 3.0 - pi^2/4.0;

  for n=1:100
    df = 3.0*l/4.0 .* (4*n+5) ./ ((2*n+3)*(n+1))^2 ./ ((2*n+3)*(n+1) - l);
    f = f + df;
    if max(df) < 1e-8; break; end
  end
end

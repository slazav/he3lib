% Yosida function integrand
function int = integrand_s(x,gap,ttc)
  c = 4;  % important power factor
  ksi = atanh(x) * c;
  Ek=sqrt(ksi.^2+gap.^2);
  phi=(cosh(Ek/(2*ttc))).^(-2) / 2 / ttc;
  int=(ksi./ttc).^2 .* phi ./ (1-x.^2) * c;
end

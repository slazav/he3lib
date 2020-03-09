% Yosida function integrand
function int = integrand_c(x,gap,dgap2,ttc)
  c = 4;  % important power factor
  ksi = atanh(x) * c;
  Ek=sqrt(ksi.^2+gap.^2);
  phi=(cosh(Ek/(2*ttc))).^(-2) / 2 / ttc;
  int=(Ek.^2./ttc - 0.5*dgap2)./ttc .* phi ./ (1-x.^2) * c;
end

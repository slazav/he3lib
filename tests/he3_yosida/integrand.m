% Yosida function integrand
function int = integrand(x,gap,TTc,n)
  c = 2;  % important power factor
  ksi = atanh(x) * c;
  Ek=sqrt(ksi.^2+gap.^2);
  phi=(cosh(Ek/(2*TTc))).^(-2) / 2 / TTc;
  int=(ksi./Ek).^n .* phi ./ (1-x.^2) * c;
end

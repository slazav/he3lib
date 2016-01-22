% Solve a cubic equation A3 x^3 + A2 x^2 + A1*x + A0 = 0
% in real numbers. Output nomber of real roots and
% sorted roots x1,x2,x3 (NaN for non-existing root).
%
% the method used:
% Calculation of Densities from Cubic Equations of State: Revisited
% Ulrich K. Deiters*, Ricardo Macias-Salina
% Ind. Eng. Chem. Res., 2014, 53 (6), pp 2529�2536
%
% I want to use it in fortran, so I use this specific programming style.
%
% Here scaling of coefficients is missing (it is needed for
% large/small coefficients). See the paper or the fortran code
% in ../math.f how to add it.

function [n, x1,x2,x3] = solve_cubic(A3,A2,A1,A0)

  if A3==0
    [n,x1,x2] = solve_quadr(A0,A1,A2);
    x3=NaN;
    return;
  end

  % normalize
  a0 = A0/A3;
  a1 = A1/A3;
  a2 = A2/A3;

  if a0==0
    x1 = 0;
    [n, x2, x3] = solve_quadr(a1,a2,1);
    if n>0 && x1>x2; t=x2; x2=x1; x1=t; end % swap x1 and x2
    if n>0 && x2>x3; t=x3; x3=x2; x2=t; end % swap x2 and x3
    n=n+1;
    return;
  end

  % inflection point
  xi = -a2/3;
  yi = xi^3 + a2*xi^2 + a1*xi + a0;
  if yi==0
    x2 = xi;
    c1 = x2+a2;
    c0 = c1*x2 + a1;
    [n, x1, x3] = solve_quadr(c0,c1,1);
    return;
  end

  D = a2^2-3*a1;
  if D==0
    x1=xi-yi^(1/3.0);
    x2=NaN;
    x3=NaN;
    n=1;
    return;
  end

  x1=xi;
  if D>0; x1=xi-yi/abs(yi) *2/3*sqrt(D); end

  while 1;
    y  = x1^3 + a2*x1^2 + a1*x1 + a0;
    yp = 3*x1^2 + 2*a2*x1 + a1;
    ypp = 6*x1 + 2*a2;
    dx = y*yp/(yp^2-0.5*y*ypp);
    if abs(dx/x1)<1e-10; break; end
    x1 = x1-dx;
  end

  if D>0;
    c1 = x1+a2;
    c0 = c1*x1 + a1;
    [n x2,x3] = solve_quadr(c0,c1,1);
    if n>0 && x1>x2; t=x2; x2=x1; x1=t; end % swap x1 and x2
    if n>0 && x2>x3; t=x3; x3=x2; x2=t; end % swap x2 and x3
    n=n+1;
    return;
  end
  n=1;
  x2=NaN;
  x3=NaN;
  return;
end

% solve quadratic equation A2 x^2 + A1*x + A0 = 0
function [n, x1, x2] = solve_quadr(A0,A1,A2)
  D = A1^2 - 4*A2*A0;
  if (D<0)
    n=0;
    x1=NaN;
    x2=NaN;
  else
    n=2;
    x1 = (-A1-sqrt(D))/(2*A2);
    x2 = (-A1+sqrt(D))/(2*A2);
  end
end


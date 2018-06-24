function plot_spec()


% solve quadratic equation A2 x^2 + A1*x + A0 = 0
function [x1, x2] = solve_quadr(A0,A1,A2)
  D = A1^2 - 4*A2*A0;
  x1 = (-A1-sqrt(D))/(2*A2);
  x2 = (-A1+sqrt(D))/(2*A2);
end

#  Equation:
#
#  - s2*x^6
#  + (1+b)*x^4
#  - 1/sqrt(15) (3/8*(sqrt(15)*d - s*h) + b*(sqrt(15)*d + 3*s*h))*x^2
#  + b/10*(sqrt(15)*d - s*h) = 0
#
#  x = fp/f0
#  s1 - direction of n
#  s2 - sign of theta
#  s = s1*s2
#  d = (f-f0)/f0
#  h = Hr/H0
#  b = (fB/f0)^2

#  - s2*x^6
#  + (1+b)*x^4
#  - 1/sqrt(15) (- 3/8*s*h) + b*3*s*h)*x^2
#  - b/10 (s*h)^2 = 0

s1 = +1;
s2 = +1;

fB = 100e3;
f0 = 1120e3;
h = 0.01;
dd=0:1e-5:0.02;

b = (fB/f0)^2;
s = s1*s2;

for i=1:length(dd)
  d=dd(i);
  s=+1;
  a2 = 1+b;
  a1 = -1/sqrt(15)*(3/8.0*(sqrt(15)*d - s*h) + b*(sqrt(15)*d + 3*s*h));
  a0 = b/10*(sqrt(15)*d - s*h)*s*h;
#  [n(i) x1(i) x2(i) x3(i)] = solve_cubic(a3,a2,a1,a0);

  [x1a(i) x2a(i)] = solve_quadr(a0,a1,a2);

  s=-1;
  a2 = 1+b;
  a1 = -1/sqrt(15)*(3/8.0*(sqrt(15)*d - s*h) + b*(sqrt(15)*d + 3*s*h));
  a0 = b/10*(sqrt(15)*d - s*h)*s*h;
  [x1b(i) x2b(i)] = solve_quadr(a0,a1,a2);
end

find_figure('plot_spec'); clf; hold on;

plot(dd, imag(sqrt(x1a)), 'b-')
plot(dd, imag(sqrt(x2a)), 'b-')
plot(dd, imag(sqrt(x1b)), 'c-')
plot(dd, imag(sqrt(x2b)), 'c-')

plot(dd, real(sqrt(x1a)), 'r-')
plot(dd, real(sqrt(x2a)), 'r-')

plot(dd, real(sqrt(x1b)), 'm-')
plot(dd, real(sqrt(x2b)), 'm-')



#plot(dd, x3, 'b-')
plot(dd, zeros(size(dd)), 'k-')

end


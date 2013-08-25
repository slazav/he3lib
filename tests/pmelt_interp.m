#!/usr/bin/octave -qf
# interpolation of melting pressure between measured ranges.

function p = pmelt1(t)
% Greywall-86
 p =  34.3380D0 ...
  - 0.19652970D-1 * (t*1D3).^(-3) ...
  + 0.61880268D-1 * (t*1D3).^(-2) ...
  - 0.78803055D-1 * (t*1D3).^(-1) ...
  + 0.13050600D0 ...
  - 0.43519381D-1 * (t*1D3) ...
  + 0.13752791D-3 * (t*1D3).^2 ...
  - 0.17180436D-6 * (t*1D3).^3 ...
  - 0.22093906D-9 * (t*1D3).^4 ...
  + 0.85450245D-12* (t*1D3).^5;
end
function dp = dpmelt1(t)
 dp = ...
  + 3 * 0.19652970D-1 * (t*1D3).^(-4) ...
  - 2 * 0.61880268D-1 * (t*1D3).^(-3) ...
  + 0.78803055D-1 * (t*1D3).^(-2) ...
  - 0.43519381D-1 ...
  + 2 * 0.13752791D-3 * (t*1D3).^1 ...
  - 3 * 0.17180436D-6 * (t*1D3).^2 ...
  - 4 * 0.22093906D-9 * (t*1D3).^3 ...
  + 5 * 0.85450245D-12* (t*1D3).^4;
 dp = dp*1000;
end

function p = pmelt1a(t)
  % plts2000
  p = ...
  - 1.3855442D-12 * t.^(-3) ...
  + 4.5557026D-9  * t.^(-2) ...
  - 6.4430869D-6  * t.^(-1) ...
  + 3.4467434D0 ...
  - 4.4176438D0 * t.^1 ...
  + 1.5417437D1 * t.^2 ...
  - 3.5789858D1 * t.^3 ...
  + 7.1499125D1 * t.^4 ...
  - 1.0414379D2 * t.^5 ...
  + 1.0518538D2 * t.^6 ...
  - 6.9443767D1 * t.^7 ...
  + 2.6833087D1 * t.^8 ...
  - 4.5875709D0 * t.^9;
  p = p * 10D0; % MPa -> bar
end
function dp = dpmelt1a(t)
  % plts2000
  dp = ...
  + 3 * 1.3855442D-12 * t.^(-4) ...
  - 2 * 4.5557026D-9  * t.^(-3) ...
  + 6.4430869D-6  * t.^(-2) ...
  - 4.4176438D0 ...
  + 2 * 1.5417437D1 * t ...
  - 3 * 3.5789858D1 * t.^2 ...
  + 4 * 7.1499125D1 * t.^3 ...
  - 5 * 1.0414379D2 * t.^4 ...
  + 6 * 1.0518538D2 * t.^5 ...
  - 7 * 6.9443767D1 * t.^6 ...
  + 8 * 2.6833087D1 * t.^7 ...
  - 9 * 4.5875709D0 * t.^8;
  dp = dp * 10D0; % MPa -> bar
end

function p = pmelt2(t)
% Osborne, Abraham, Weinstock, 1951
  p = (26.8D0 + 13.1 * t.^2) * 1.01325;
end
function dp = dpmelt2(t)
  dp = 2 * 13.1 * t * 1.01325;
end

function p = pmelt3(t)
% Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955))
 p = (25.16 + 20.08201 * t.^1.517083) * 0.980665;
end
function dp = dpmelt3(t)
 dp = 1.517083 * 20.08201 * t.^(1.517083-1) * 0.980665;
end

t1 = 0.25;
t2 = 0.50;
t3 = 1.50;
t4 = 2.00;
tm = 0.31524;
pm = 29.3113;

f1 = pmelt1(t1)
f1a = pmelt1a(t1)
f2 = pmelt2(t2)
f3 = pmelt2(t3)
f4 = pmelt3(t4)

df1 = dpmelt1(t1)
df1a = dpmelt1a(t1)
df2 = dpmelt2(t2)
df3 = dpmelt2(t3)
df4 = dpmelt3(t4)

A = [ t1^5 t1^4 t1^3 t1^2 t1 1
      t2^5 t2^4 t2^3 t2^2 t2 1
      tm^5 tm^4 tm^3 tm^2 tm 1
      5*t1^4 4*t1^3 3*t1^2 2*t1 1 0
      5*t2^4 4*t2^3 3*t2^2 2*t2 1 0
      5*tm^4 4*tm^3 3*tm^2 2*tm 1 0];
B = [f1; f2; pm; df1; df2; 0];
Ba = [f1a; f2; pm; df1a; df2; 0];

fprintf('>> %f %f\n', t1, t2);
X = inv(A) * B;
fprintf('%15f\n', X);

fprintf('>> plts %f %f\n', t1, t2);

X = inv(A) * Ba;
fprintf('%15f\n', X);

A = [ t3^3 t3^2 t3 1
      t4^3 t4^2 t4 1
      3*t3^2 2*t3 1 0
      3*t4^2 2*t4 1 0];
B = [f3; f4; df3; df4];

fprintf('>> %f %f\n', t3, t4);
X = inv(A) * B;
fprintf('%15f\n', X);


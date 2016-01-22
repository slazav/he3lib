function plot_phase()
% b2 phase magnetization and susceptibility

  figure(1); clf; hold on;

  pc=34.338;
  P=(0.0:0.2:33);
  p=P/pc;
  p2 = p.^2;
  p3 = p.^3;

  A(1) =   3391; A(2) =   21500; A(3) =  -8490; A(4) =       0; A(5) = 2.098 ; % Bc
  B(1) =   1.41; B(2) =       0; B(3) =      0; B(4) =       0; B(5) =      0; % f3
  C(1) =  -0.29; C(2) =   -0.41; C(3) =      0; C(4) =       0; C(5) =      0; % f4
  D(1) = 0.9294; D(2) =   6.659; D(3) = 0.3406; D(4) = -0.5012; D(5) = 1.9844; % Tc
  E(1) = 1.5745; E(2) = -1.1222; E(3) = 0.3242; E(4) =       0; E(5) =      0; % Tab/Tc
  F(1) =  0.616; F(2) = -1.174 ; F(3) = 0.301 ; F(4) =       0; F(5) =      0; % 1/g
  G(1) =  0.303; G(2) =  0.717 ; G(3) = 0.112 ; G(4) =       0; G(5) = 3.611 ; % 1+Foa

  Bc  = (A(1) + A(2)*p + A(3)*p2 + A(4)*p3) ./ (1 + A(5)*p );
  f3  = (B(1) + B(2)*p + B(3)*p2 + B(4)*p3) ./ (1 + B(5)*p );
  f4  = (C(1) + C(2)*p + C(3)*p2 + C(4)*p3) ./ (1 + C(5)*p );
  Tc  = (D(1) + D(2)*p + D(3)*p2 + D(4)*p3) ./ (1 + D(5)*p );
  Tab = (E(1) + E(2)*p + E(3)*p2 + E(4)*p3) ./ (1 + E(5)*p );
  gr  = (F(1) + F(2)*p + F(3)*p2 + F(4)*p3) ./ (1 + F(5)*p );
  Foa = (G(1) + G(2)*p + G(3)*p2 + G(4)*p3) ./ (1 + G(5)*p );

  Foa = he3_f0a(P)+1;

  Bo  = 1.97*Tc.*Foa*10000;
  BB  = gr/4.*(Bo./Bc).^2;

  % original f5
  f5  = (1-f3.*Tab.^6-f4.*Tab.^8-(1-f3-f4).*Tab.^2 ...
       +(1+2*f3+3*f4).*(Tab.^4-Tab.^2))./ ...
        (BB .* (Tab.^4-Tab.^2))-1;

  %% fit f5
  %a = fit(p', f5', '(a + b*x + c*x^2 + d*x^3 + e*x^4)/(1+f*x)', 'startpoint', [1 1 1 1 1 1])
  %fprintf('%f %f %f %f %f %f\n', a.a,a.b,a.c,a.d,a.e,a.f);
  %plot(P, a(p), 'm-')

  %my own f5
  f5 = (0.041870 + 5.417531*p -10.044312*p2 + 7.639438*p3 - 2.379517*p2.*p2)./(1-0.537422*p);
  %plot(P, f5, 'r-')
  %plot(P, f5a, 'g-')

  f2  = BB.*(1+f5)-(1+2*f3+3*f4);
  f1  = 1+f5-f2-f3-f4;
  fn  = @(x) (f1*x + f2*x^2 + f3*x^3 + f4*x^4)./(1 + f5*x);

  %plot(P, Tc, 'r-')
  %plot(P, Tc.*Tab, 'b-')
  %return

  for ttc=0:0.01:1
    H=he3_b2hcr(ttc,p);
    x2 = fn(ttc^2);
    H = sqrt(1-x2.^2).*Bc;
    plot(P, H, 'r-')
    plot(P, he3_b2hcr(ttc,P), 'g-')
  end
%  plot(P, Bc, 'b-')
end

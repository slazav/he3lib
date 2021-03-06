function plot_b2rhon()
% This is calculation of B-phase normal fluid density

  figure(1);
  clf; hold on;

  p=0;


  if (0) % plot integrand
    ttc=0.05;
    gap=he3_gap(ttc,p);
    x=0.001:0.01:0.999;
    ct=0:0.01:1;
    [xx,yy]=meshgrid(x,ct);
    int = integrand(xx,yy, gap*1.1,gap*0.9,ttc, 1);
    surface(xx,yy,int,'EdgeColor','none');
  end

  if (0)
      ttc=0.01:0.03:1;
      for i=1:length(ttc);
        H = he3_b2hcr(ttc(i),p);
        [rpar0(i) rper0(i)] = rho(ttc(i),p,0);
        [rpar(i) rper(i)]   = rho(ttc(i),p,H);
      end
      plot(ttc, rpar, 'r.-')
      plot(ttc, rper, 'b.-')

      plot(ttc, rpar0, 'm.-')
      plot(ttc, rper0, 'c.-')
      plot(ttc, he3_rho_nb(ttc,p), 'k.-')
  end

  % test the library function
  if (0)
      ttc=0.01:0.03:1;
      H = he3_b2hcr(ttc,p);
      rpar0 = he3_b2rho_npar(ttc,p,0);
      rper0 = he3_b2rho_nper(ttc,p,0);
      rpar = he3_b2rho_npar(ttc,p,H);
      rper = he3_b2rho_nper(ttc,p,H);

      plot(ttc, rpar, 'r.-')
      plot(ttc, rper, 'b.-')
      plot(ttc, rpar0, 'm.-')
      plot(ttc, rper0, 'c.-')
      plot(ttc, he3_rho_nb(ttc,p), 'k.-')
      title(['Normal fluid density at the A-B boundary, P = ', num2str(p)]);
      xlabel('T/Tc')
      ylabel('\rho/\rho_n')
      legend('\rho_{par}', '\rho_{perp}', '\rho_n')
  end

  % test the library function -2
  if (1)
      ttc=0.07:0.03:1;
      H = he3_b2hcr(ttc,p);
      rpar0 = he3_b2rho_npar(ttc,p,0);
      rper0 = he3_b2rho_nper(ttc,p,0);
      rpar = he3_b2rho_npar(ttc,p,H);
      rper = he3_b2rho_nper(ttc,p,H);

      plot(ttc, rpar./rpar0, 'r.-')
      plot(ttc, rper./rper0, 'b.-')
      plot([0 1], [1 1], 'k')
      ylim([0 3])
      title(['Normal fluid density at the A-B boundary, P = ', num2str(p)]);
      xlabel('T/Tc')
      ylabel('\rho/\rho(H=0)')
      legend('\rho_{par}', '\rho_{perp}', '\rho_n')
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function int = integrand(x,ct,gap1,gap2,ttc, n)

  c = 2;  % important power factor
  xi = atanh(x) * c;
  dxi = 1./(1-x.^2) * c;
  gg = gap1^2*(1-ct.^2) + gap2^2*ct.^2;
  Ek = sqrt(xi.^2+gg);

  % df/dEk
  phi=-(cosh(Ek/(2*ttc))).^(-2)/4/ttc;

  if n==1
    int= - 2*phi .* (ct.^2) .* dxi;
  else
    int= - 2*phi .* (1-ct.^2)/2 .* dxi;
  end
end

function [rpar rper] = rho(ttc,p,H)
  gap1 = he3_b2gap1(ttc,p,H);
  gap2 = he3_b2gap2(ttc,p,H);
  if any([~isfinite(gap1) ~isfinite(gap2)])
    rpar=NaN; rper=NaN; return;
  end
  ff1 = @(x,y) integrand(x,y,gap1,gap2, ttc, 1);
  ff2 = @(x,y) integrand(x,y,gap1,gap2, ttc, 2);
  rpar = quad2d(ff1,0,1,0,1)*3;
  rper = quad2d(ff2,0,1,0,1)*3;
  % fermi-liquid corrections:
  f1s = he3_f1s(p);
  rpar = rpar * (1+f1s/3)/(1+f1s/3*rpar);
  rper = rper * (1+f1s/3)/(1+f1s/3*rper);
end

function process_data()
  pkg load optim
  find_figure('visc'); clf; hold on;

  global dofit
  dofit = 0;  % 0 - no fit / 1 - 2-par fit / 2 - 1-par fit (alpha only)

  k = 0.893 # Alvesalo -> Greywall temperature scale

  p=[]; EA=[]; EG=[];
  [p(end+1), EG(end+1), EA(end+1)] = process_file('../1996_nakagawa_visc/fig1_00.txt', 0.0, 1, 'ro');
  [p(end+1), EG(end+1), EA(end+1)] = process_file('../1996_nakagawa_visc/fig1_05.txt', 5.0, 1, 'go');
  [p(end+1), EG(end+1), EA(end+1)] = process_file('../1996_nakagawa_visc/fig1_21.txt', 21.0, 1,  'bo');
  [p(end+1), EG(end+1), EA(end+1)] = process_file('../1996_nakagawa_visc/fig1_29.txt', 29.0, 1, 'mo');

  [p(end+1), EG(end+1), EA(end+1)] = process_file('fig4_02.10.txt', 02.10, k, 'b*');
  [p(end+1), EG(end+1), EA(end+1)] = process_file('fig4_04.65.txt', 04.65, k, 'm*');
  [p(end+1), EG(end+1), EA(end+1)] = process_file('fig4_09.89.txt', 09.89, k, 'c*');
  [p(end+1), EG(end+1), EA(end+1)] = process_file('fig4_19.89.txt', 19.89, k, 'r*');
  [p(end+1), EG(end+1), EA(end+1)] = process_file('fig4_29.34.txt', 29.34, k, 'g*');

  % black curves - fits by CHH
  if (0)
    PCHH = [0.1 1.28 2.10 4.65 9.89 19.89 29.34];
    ACHH = [0.373 0.380 0.382 0.424 0.505 0.603 0.710];
    BCHH = [0.04 0.06 0.09 0.22 0.38 1.13 1.45];
    for i=1:length(PCHH)
      tc = he3_tc(PCHH(i));
      tt = linspace(tc, 9,30);
      vi = 1./(ACHH(i)+BCHH(i)./(tt/k).^2) * k^2;
      plot(tt, vi, 'k-');
    end
  end

  fa = @(pp,x)  pp(1) + pp(2)./(x + pp(3));
  px = 0:30;
  pa = nonlin_curvefit(fa, [1;1;1], p, EA);

  find_figure('emery - alpha'); clf; hold on;
  plot(p, EA, 'r*');
  plot(px, fa(pa,px))

  pg = mean(EG);

  find_figure('emery - G'); clf; hold on;
  plot(p, EG, 'r*');
  plot([0 30], [pg pg]);

  printf('pg=%e;\n', pg)
  printf('pa=[%e %e %e];\n', pa)

end

function [p, EG,EA] = process_file(fname, p, k, c)
  global dofit

  ## read data file
  ff = fopen(fname);
  data = textscan(ff, '%f %f\n', 'CommentStyle', '#');
  T = data{1}; # T
  E = data{2}; # eta*T^2

  # temperature correction
  Tc = he3_tc(p);
  T = T*k;
  E = E*k^2;
  ii = find(T>Tc);
  T = T(ii); E = E(ii);

  # viscosity from he3lib
  E0 = visc0(T*1e-3,p).*T.^2*1e-6;


  # Emery function
  ff = @(pp,ttc) 1 - pp(1)*(1 - sqrt(ttc-1).*atan(pp(2)./sqrt(ttc-1))/pp(2));

  # initial values - fit result
  pg=2.318858e-01;
  pa=[-4.067317e-02 4.841729e+01 6.0755];
  pp= [pg; pa(1) + pa(2)/(p + pa(3))];

  # fit, or use previously fitted model:
  if (dofit == 1)
    pp = nonlin_curvefit(ff, pp, T/Tc, E./E0);

  elseif (dofit == 2)
    f1 = @(pp1,x) ff([pp(1) pp1], x);
    pp(2) = nonlin_curvefit(f1, pp(2), T/Tc, E./E0);

  end

  tt = linspace(Tc,9,50);
  EE0 = visc0(tt*1e-3,p).*tt.^2*1e-6;
  EF = EE0.*ff(pp,tt/Tc);

  plot (T,E, c)
  plot (tt,EF, [c(1) '-'])
  plot (tt,EE0, [c(1) '--'])
  text (tt(end)+0.1, EF(end), sprintf('p = %.2f bar', p))
  printf('%f %f %f\n', p, pp)

  EG = pp(1);
  EA = pp(2);
end


# Dyugaev model for viscosity
function E = visc0(t,p)
  Te = 0.064795 + p*2.949606e-5 + 5.351244./(p+16.87556);
  E0 = 22.3125  + p*0.375 - 36.09375./(p+7.5);
  E = E0.*((Te./t).^2 + 1.41*Te./t + 1);
end


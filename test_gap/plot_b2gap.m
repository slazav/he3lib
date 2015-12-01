function plot_b2gap1()
% This is calculation of B-phase gap distortion and spin polarization,
% based on Ashida and Nagai paper (Progr.Theor.Phys. 74 949 (1985)).
%
% The free energy is written as a function af three parameters:
%  * gap_perp = gap1 = gap*(1+A)
%  * gap_parallel = gap2 = gap*(1-B)
%  * effective field we
% This program calculates the minimum of the energy


  figure(1);
  clf; hold on;
  p=0.523;     % tc=1mK, as in the paper
  f0a = he3_f0a(p);
  f0a = -0.75; % as in the paper

  ttc=0.1;
  gap=he3_gap(ttc,p);

  xx=-1:0.1:1;
  yy=-1:0.1:1;
  zz = 0.2;
  w0 = 0.1;

%  x=0.001:0.001:0.999;
%  ct=0:0.001:1;
%  [xx,yy]=meshgrid(x,ct);
%  int = integrand(xx,yy, gap,gap*1.01,gap*0.99,ttc, 0.05, 1);
%  surface(xx,yy,int,'EdgeColor','none');

%  xx=0:0.01:3;
%  w0=0.1; A=0.0052; B=-0.0110;
%  for i=1:length(xx)
%    int1(i) = -zerofunc1(ttc, gap, f0a, w0, xx(i), A, B);
%  end
%  tc = he3_tc(p)*1e-3*const_kb; % erg units
%  xxh = xx*gap*tc/he3_gyro/const_hbar; % G
%  int1 = int1*tc *he3_2n0(p)/2; % dimensionless
%
%  mag = int1*he3_gyro*const_hbar;
%
%  chi1=he3_chi_n(p)*he3_chi_b(p,ttc);
%  mag1=chi1*xxh;

%  plot(xxh, mag,'b-');
%  plot(xxh, mag1,'r-');


  w0=0.0:0.01:0.4;
  figure(2);
  clf; hold on;

  res=calc(w0,f0a,p,ttc);
  plot(res.w0,res.A, 'r-')
  plot(res.w0,res.B, 'b-')

  % table for Fischer's thesis for comparison
  f_h  = 0.05:0.01:0.33;
  f_D1 = [13 20 26 34 43 52 63 75 88 102 118 ...
          134 152 171 191 213 237 262 288 314 ...
          345 377 412 450 490 534 584 640 705]*1e-4 + 1;
  f_D2 = [997 996 995 993 991 989 987 985 982 ...
          979 976 972 968 964 959 954 949 943 ...
          936 930 922 914 905 894 883 870 855 837 814]*1e-3;
  plot(f_h,f_D1-1, 'ro')
  plot(f_h,f_D2-1, 'bo')
end

function res = calc(w0,f0a,p,ttc);
 he3_gyro=2.037800e+04;
  gap=he3_gap(ttc,p);
  tc = he3_tc(p)*1e-3*const_kb; % erg units
  w0G = w0*gap*tc/he3_gyro/const_hbar; % G
  for i=1:length(w0)
    % initial values:
    we(i)=(1-4*f0a*he3_chi_b(ttc,p))*w0(i);
    A(i)=0.1; B(i)=-0.1;
    [we(i) A(i) B(i)] = find_min(w0(i),f0a,p,ttc);
  end
  res.w0=w0;
  res.w0G=w0G;
  res.A=A;
  res.B=B;
end

function [we A B] = find_min(w0,f0a,p,ttc);
  gap=he3_gap(ttc,p);
  ff = @(x) zerofunc(ttc, gap, f0a, w0, x(1),x(2),x(3));

  x=[2*w0,0,0];
  f1 = ff(x);
  sl = 3; % slope
  for i=1:100
    dx = sl*f1;
    x  = x+dx;
    f2 = ff(x);

    if sum(f2.^2)<1e-12; break; end
    if abs(x(2))>1 || abs(x(3))>1 || sl<0;
      x=[NaN NaN NaN]; break; end
    sl = dx/(f1-f2);
    f1=f2;
  end
  fprintf('> %f %f %f %d %f\n', x, i, sl);
  we=x(1);
  A=x(2);
  B=x(3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% integrands for dF/dgap1^2, dF/dgap2^2, dF/we
function int = integrand(x,ct,gap0,gap1,gap2,TTc,we, n)
  xi = 1./x-1;

  gaps = gap1^2*(1-ct.^2) + gap2^2*ct.^2;
  Ez = sqrt(xi.^2 + gap2^2*ct.^2);
  EB = sqrt(xi.^2 + gap0^2);
  Esp = @(s) sqrt(xi.^2 + gaps + we^2/4 - s*we*Ez);

  % theta/2E
  Ap = -1/4*tanh(Esp(1)/(2*TTc))./Esp(1);
  Am = -1/4*tanh(Esp(-1)/(2*TTc))./Esp(-1);
  AB = -1/4*tanh(EB/(2*TTc))./EB;

  if n==1;
    int  = (Ap.*(we/2-Ez) + Am.*(we/2+Ez))/2;
  elseif n==2;
    int  = (Ap+Am-2*AB).*(1-ct.^2);
  elseif n==3;
    int  = (Ap.*(1-we./(2*Ez)) + Am.*(1+we./(2*Ez))-2*AB).*ct.^2;
  end

  int = - int ./ x.^2 * 2;
  int(find(isnan(int)))=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions for minimization
%% ttc - T/Tc
%% gap - gap/Tc
%% f0a -- Fermi-liquid parameter
%% A,B -- gap distortions
%% w0,we -- external and effective fields, in gap units
function v = zerofunc(ttc, gap, f0a, w0, we, A,B)
  we=we*gap; w0=w0*gap;
  gap1 = gap*(1+A);
  gap2 = gap*(1+B);
  f1 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,1);
  f2 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,2);
  f3 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,3);
  v1 = quad2d(f1,0,1,0,1) + (we-w0)/(4*f0a);
  v2 = quad2d(f2,0,1,0,1);
  v3 = quad2d(f3,0,1,0,1);
  %  fprintf('%f %f %f %f -- %f %f %f\n', w0, A, B, we, v1,v2,v3)
  v = [v1 v2 v3];
end

% Integration. I have no quad2 function in octave.
% Integrand is good, no need to use any smart integration here.
function v=quad2d(f, x1,x2,y1,y2)
  n=30;
  dx=(x2-x1)/(n-1);
  dy=(y2-y1)/(n-1);
  [xx,yy] = meshgrid(x1:dx:x2,y1:dy:y2);
  v=sum(sum(f(xx,yy)))*dx*dy;
end

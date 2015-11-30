function plot_g1_int()
% This is calculation of B-phase gap distortion and spin polarization,
% based on Ashida and Nagai paper (Progr.Theor.Phys. 74 949 (1985)).
%
% The free energy is written as a function af three parameters:
%  * gap_perp = gap1 = gap*(1+A)
%  * gap_parallel = gap2 = gap*(1-B)
%  * effective field we
% This program calculates the minimum of the energy

  find_figure('b2gap');
  clf; hold on;
  p=0.523;     % tc=1mK, as in the paper
  %f0a = he3_f0a(p);
  f0a = -0.75; % as in the paper

  ttc=0.1;
  gap=he3_gap(ttc,p);


%  x=0.001:0.001:0.999;
%  ct=0:0.001:1;
%  [xx,yy]=meshgrid(x,ct);
%  int = integrand(xx,yy, gap,gap*1.01,gap*0.99,ttc, 0.05, 1);
%  surface(xx,yy,int,'EdgeColor','none');

%  % simple test -- find magnetization for a fixed gap distortion A,B
%  % as a function of effective field
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


  w0=0.03:0.03:0.33;
  find_figure('b2gap 1');
  clf; hold on;

  res=calc(w0,f0a,p,ttc);
  plot(res.w0,res.A, 'r.-')
  plot(res.w0,res.B, 'b.-')

  % table for Fischer's thesis for comparison
  f_h  = 0.05:0.01:0.33;
  f_D1 = [13 20 26 34 43 52 63 75 88 102 118 ...
          134 152 171 191 213 237 262 288 314 ...
          345 377 412 450 490 534 584 640 705]*1e-4 + 1;
  f_D2 = [997 996 995 993 991 989 987 985 982 ...
          979 976 972 968 964 959 954 949 943 ...
          936 930 922 914 905 894 883 870 855 837 814]*1e-3;
  plot(f_h,f_D1-1, 'r--')
  plot(f_h,f_D2-1, 'b--')
end

function res = calc(w0,f0a,p,ttc);
  gap=he3_gap(ttc,p);
  tc = he3_tc(p)*1e-3*const_kb; % erg units
  w0G = w0*gap*tc/he3_gyro/const_hbar; % G

  for i=1:length(w0)
    % initial values:
    we(i)=(1-4*f0a*he3_chi_b(ttc,p))*w0(i);
    A(i)=0.1; B(i)=-0.1;
    vt=0;
    if (vt==0) % minimize all three functions
      func=@(pp) zerofunc(ttc, gap, f0a, w0(i), pp(1),pp(2),pp(3));
      pp = fminsearch(func, [we(i), A(i), B(i)]);
      we(i)=pp(1);
      A(i)=pp(2);
      B(i)=pp(3);
    elseif (vt==1) % calculate spin polarization from known susceptibility
       we(i)= 2.2*w0(i);
      func=@(pp) zerofunc2(ttc, gap, f0a, w0(i), we(i), pp(1),pp(2));
      pp = fminsearch(func, [A(i), B(i)]);
      A(i)=pp(1);
      B(i)=pp(2);
    end
    fprintf('%f %f  %f %f  %f\n', w0(i), w0G(i), A(i), B(i), we(i)/w0(i))
  end
  res.w0=w0;
  res.w0G=w0G;
  res.A=A;
  res.B=B;
end

function int = integrand(x,ct,gap0,gap1,gap2,TTc,we, n)
  %c = 2;  % important power factor
  %xi = atanh(x) * c;
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
    int  = -(Ap.*(we/2-Ez) + Am.*(we/2+Ez));
    int0 = 0;
  elseif n==2;
    int  = (Ap+Am).*(1-ct.^2);
    int0 = 2*AB.*(1-ct.^2);
  elseif n==3;
    int  = (Ap.*(1-we./(2*Ez)) + Am.*(1+we./(2*Ez))).*ct.^2;
    int0  = 2*AB.*ct.^2;
  end

  int = - (int-int0) ./ x.^2 * 2;
  int(find(isnan(int)))=0;
end

%% gap in tc units, w0,we in gap units!
function v = zerofunc(ttc, gap, f0a, w0, we, A,B)
  we=we*gap; w0=w0*gap;
  gap1 = gap*(1+A);
  gap2 = gap*(1+B);
  f1 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,1);
  f2 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,2);
  f3 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,3);
  v1 = quad2d(f1,0,1,0,1)/2 - (we-w0)/(4*f0a);
  v2 = quad2d(f2,0,1,0,1);
  v3 = quad2d(f3,0,1,0,1);
  %  fprintf('%f %f %f %f -- %f %f %f\n', w0, A, B, we, v1,v2,v3)
  v = v1^2 + v2^2 + v3^2;
end

%% gap, w0,we in tc units!
function v = zerofunc2(ttc, gap, f0a, w0, we, A,B)
  we=we*gap; w0=w0*gap;
  gap1 = gap*(1+A);
  gap2 = gap*(1+B);
  f2 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,2);
  f3 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,3);
  v2 = quad2d(f2,0,1,0,1);
  v3 = quad2d(f3,0,1,0,1);
  %  fprintf('%f %f %f %f -- %f %f %f\n', w0, A, B, we, v1,v2,v3)
  v=v2^2+v3^2;
end

function v = zerofunc1(ttc, gap, f0a, w0, we, A,B)
  we=we*gap; w0=w0*gap;
  gap1 = gap*(1+A);
  gap2 = gap*(1+B);
  f1 = @(x,y) integrand(x,y, gap,gap1,gap2,ttc,we,1);
  v1 = quad2d(f1,0,1,0,1);
  v=v1;
end


function plot_g1_int()
% It should not go to infinity at any T/Tc.
% This can be controlled by changing c in
% x = tanh(xi/c(ttc))
% c = 2 is a rather good value

  find_figure('g1 integrand');
  clf; hold on;
%  addpath ../matlab

  p=0;
  ttc=0.1;
  gap=he3_gap(ttc,p);
  f0a = he3_f0a(p);

%  x=0.001:0.001:0.999;
%  ct=0:0.001:1;
%  [xx,yy]=meshgrid(x,ct);
%  int = integrand(xx,yy, gap,gap*1.01,gap*0.99,ttc, 0.05, 1);
%  surface(xx,yy,int,'EdgeColor','none');

%  xx=0.:0.01:3;
%  w0=0.2;
%  for i=1:length(xx)
%    int1(i) = zerofunc(ttc, gap, f0a, w0*gap, xx(i)*gap, 1.0213, 0.954);
%    int2(i) = zerofunc(ttc, gap, f0a, w0*gap, xx(i)*gap, 1.001339, 0.99266);
%  end
%  plot(xx,int1,'r-');
%  plot(xx,int2,'b-');

  w0=0.03:0.03:0.33;
  for i=1:length(w0)
    % initial values:
    we(i)=(1-4*f0a*he3_chi_b(ttc,p))*w0(i);
    A(i)=1.01; B(i)=0.99;
    vt=0;
    if (vt==0) % minimize all three functions
      func=@(pp) zerofunc(ttc, gap, f0a, w0(i), pp(1),pp(2),pp(3));
      pp = fminsearch(func, [we(i), A(i), B(i)]);
      we(i)=pp(1);
      A(i)=pp(2);
      B(i)=pp(3);
    elseif (vt==1) % calculate spin polarization from known susceptibility
      func=@(pp) zerofunc2(ttc, gap, f0a, w0(i), we(i), pp(1),pp(2));
      pp = fminsearch(func, [A(i), B(i)]);
      A(i)=pp(1);
      B(i)=pp(2);
    else % "smart" iterative method
      for j=1:4
        func=@(pp) zerofunc2(ttc, gap, f0a, w0(i), we(i), pp(1),pp(2));
        pp = fminsearch(func, [A(i), B(i)]);
        A(i)=pp(1);  B(i)=pp(2);

        func=@(pp) zerofunc1(ttc, gap, f0a, w0(i), pp, A(i),B(i));
        pp = fminsearch(func, we(i));
        we(i) = pp;
      end
    end
    fprintf('>> %f %f %f %f\n', w0(i), A(i), B(i), we(i))
  end

  find_figure('gap distortion');
  clf; hold on;
  plot(w0,(A-1), 'r.-')
  plot(w0,(B-1), 'b.-')
%  print -deps -color plot_yosida_int.eps

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
    int  = Ap.*(we/2-Ez) + Am.*(we/2+Ez);
    int0 = 0;
  elseif n==2;
    int  = (Ap+Am).*(1-ct.^2);
    int0 = 2*AB.*(1-ct.^2);
  elseif n==3;
    int  = (Ap.*(1-we./(2*Ez)) + Am.*(1+we./(2*Ez))).*ct.^2;
    int0  = 2*AB.*ct.^2;
  end

  int = - (int-int0) ./ x.^2 * 2;
end

%% gap in tc units, w0,we in gap units!
function v = zerofunc(ttc, gap, f0a, w0, we, A,B)
  we=we*gap; w0=w0*gap;
  f1 = @(x,y) integrand(x,y, gap,gap*A,gap*B,ttc,we,1);
  f2 = @(x,y) integrand(x,y, gap,gap*A,gap*B,ttc,we,2);
  f3 = @(x,y) integrand(x,y, gap,gap*A,gap*B,ttc,we,3);
  v1 = quad2d(f1,0,1,0,1) + 2*(we - w0)/(f0a/2)^2 - (w0-2*we)/(f0a/2);
  v2 = quad2d(f2,0,1,0,1);
  v3 = quad2d(f3,0,1,0,1);
  %  fprintf('%f %f %f %f -- %f %f %f\n', w0, A, B, we, v1,v2,v3)
  v = v1^2 + v2^2 + v3^2;
end

%% gap, w0,we in tc units!
function v = zerofunc2(ttc, gap, f0a, w0, we, A,B)
  we=we*gap; w0=w0*gap;
  f2 = @(x,y) integrand(x,y, gap,gap*A,gap*B,ttc,we,2);
  f3 = @(x,y) integrand(x,y, gap,gap*A,gap*B,ttc,we,3);
  v2 = quad2d(f2,0,1,0,1);
  v3 = quad2d(f3,0,1,0,1);
  %  fprintf('%f %f %f %f -- %f %f %f\n', w0, A, B, we, v1,v2,v3)
  v=v2^2+v3^2;
end

function v = zerofunc1(ttc, gap, f0a, w0, we, A,B)
  we=we*gap; w0=w0*gap;
  f1 = @(x,y) integrand(x,y, gap,gap*A,gap*B,ttc,we,1);
  v1 = quad2d(f1,0,1,0,1)*4*f0a + we - w0;
  %  fprintf('%f %f %f %f -- %f %f %f\n', w0, A, B, we, v1,v2,v3)
  v=v1^2;
end


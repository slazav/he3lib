function plot_kspec()
  % original program for plotting B-phase spin-wave spectra
  % (full equation for a uniform texture)
  % use two solvers: strightforward non-linear one, and cubic equation solver
  % (which was transferred then into fortran library)

  figure(1); clf; hold on;

  bn  = 40; % beta_n, deg
  an  =  0; % alpha_n, deg
  bk  =  90; % beta_k, deg
  ak  = 0;  % alpha_k, deg
  wL  = 1.5;  % wL/wB
  c12 = 2;  % gradient term, K
  c22 = 1;  % gradient term, K'
  km  = 2; % maximal k value

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  n  = [sind(bn)*cosd(an); sind(bn)*sind(an); cosd(bn)];
  kn = [sind(bk)*cosd(ak); sind(bk)*sind(ak); cosd(bk)];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function R=rmatrix(n)
    nn = kron(n,n');
    en = [    0   n(3) -n(2)
           -n(3)    0   n(1)
            n(2) -n(1)    0];
    dd=diag([1,1,1]);
    th=acos(-1/4);
    R = dd*cos(th) + nn*(1-cos(th)) - en*sin(th);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% equation for frequency in the form det(A(w,k))=0
  function y=specfunc(w,wL,k,n);

    %%% R matrix
    nn = kron(n,n');
    dd=diag([1,1,1]);
    Rk = rmatrix(n)*k;

    A = dd*w^2;
    A = A + [      0, -1i*w*wL, 0
             1i*w*wL,        0, 0
                   0,        0, 0];

    %%% dipolar and gradient torque
    A = A - nn - c12*dd*sum(k.^2) + c22*kron(Rk,Rk');
    y=real(det(A));
  end

  kv=linspace(0,km,100);

  % accurate values for k=0
  cb = n(3);
  w0(1) = 0;
  w0(2) = sqrt( (wL^2+1)/2 - sqrt( (wL^2+1)^2/4 - wL^2*cb^2 ));
  w0(3) = sqrt( (wL^2+1)/2 + sqrt( (wL^2+1)^2/4 - wL^2*cb^2 ));

  for i=1:length(kv)

    %theoretical values:
    k=kn*kv(i);
    Rk=rmatrix(n)*kn;
    Rt = (1-Rk(3)^2)/2;
    Rl = Rk(3)^2;
    w10(i)=- wL/2 + sqrt((wL/2)^2 + kv(i)^2*(c12 - c22*Rt) + 1/2*(1-cb^2));
    w20(i)= sqrt(cb^2 + kv(i)^2*(c12 - c22*Rl));
    w30(i)=  wL/2 + sqrt((wL/2)^2 + kv(i)^2*(c12 - c22*Rt) + 1/2*(1-cb^2));

    F=@(w) specfunc(w,wL,k,n);

    if 0
      % non-linear solver
      w1(i) = abs(fzero(F, w0(1)));
      w2(i) = abs(fzero(F, w0(2)));
      w3(i) = abs(fzero(F, w0(3)));

      % update initial guess for roots
      w0=[w1(i) w2(i) w3(i)];
      if i>1; w0 = 2*w0 - [w1(i-1) w2(i-1) w3(i-1)]; end
    else
      % qubic solver
      [w1a(i),w2a(i),w3a(i)] = specfunc1(c12,c22,wL,1, kv(i), ak, bk, an, bn);
      [w1b(i),w2b(i),w3b(i)] = specfunc2(c12,c22,wL,1, kv(i), ak, bk, an, bn);
    end
  end


  if (1)
    plot(0, wL,'g*', 'linewidth', 2);
    plot(0, 1 ,'g*', 'linewidth', 2);

    plot(kv, w1a,'g', 'linewidth', 2);
    plot(kv, w2a,'g', 'linewidth', 2);
    plot(kv, w3a,'g', 'linewidth', 2);

    plot(kv, w1b,'r', 'linewidth', 2);
    plot(kv, w2b,'r', 'linewidth', 2);
    plot(kv, w3b,'r', 'linewidth', 2);

    % compare specfunc1 and specfunc2
    printf('%e\n', sum(w1a-w1b).^2 + sum(w2a-w2b).^2 + sum(w3a-w3b).^2)

    plot(kv, w10,'b--');
    plot(kv, w20,'b--');
    plot(kv, w30,'b--');
    xlabel('k')
    ylabel('omega/omegaL')
  end


  % inverse spec
  if (bk!=90) return; end

  wv=linspace(1e-6,4,500);
  for i=1:length(wv)
    [k2a(i), k2b(i), k2c(i)] = specfunc_k(c12,c22, wL, 1, wv(i), an,bn);
  end

%  [wa, k2a] = fixk(wv, k2a);
%  [wb, k2b] = fixk(wv, k2b);
%  [wc, k2c] = fixk(wv, k2c);

  iia = find(k2a >=0 & k2a < km^2);
  iib = find(k2b >=0 & k2b < km^2);
  iic = find(k2c >=0 & k2c < km^2);

  plot(sqrt(k2a(iia)), wv(iia), 'g-', 'linewidth', 2)
  plot(sqrt(k2b(iib)), wv(iib), 'c-', 'linewidth', 2)
  plot(sqrt(k2c(iic)), wv(iic), 'm-', 'linewidth', 2)


end


% function sutable for fortran
%function [w1,w2,w3] = specfunc1(ttc, P, H, k, ak, bk, an, bn);
function [w1,w2,w3] = specfunc1(c1,c2,wL,wB2, kv, ak, bk, an, bn);

  % components of the n vector
  nx = sind(bn)*cosd(an);
  ny = sind(bn)*sind(an);
  nz = cosd(bn);

  % components of the k vector
  kx = kv*sind(bk)*cosd(ak);
  ky = kv*sind(bk)*sind(ak);
  kz = kv*cosd(bk);

  % rotated k vector, R_{aj} k_j
  ct=-1/4; st=sqrt(15)/4;
  kn  = kx*nx + ky*ny + kz*nz; % (k*n)
  Rkx = kx*ct + nx*kn*(1-ct) - (ky*nz-kz*ny)*st;
  Rky = ky*ct + ny*kn*(1-ct) - (kz*nx-kx*nz)*st;
  Rkz = kz*ct + nz*kn*(1-ct) - (kx*ny-ky*nx)*st;

  % he3 parameters - not now
%  cpar2 = he3_cpar(ttc, P)^2;
%  cper2 = he3_cperp(ttc, P)^2;
%  wB2 = (2*pi*he3_nu_b(ttc, P))^2;
%  wL = he3_gyro*H;

  % components of the L matrix (symmetric)
  Lxx = - c1*kv^2 + c2*Rkx*Rkx - wB2*nx*nx;
  Lyy = - c1*kv^2 + c2*Rky*Rky - wB2*ny*ny;
  Lzz = - c1*kv^2 + c2*Rkz*Rkz - wB2*nz*nz;
  Lxy =           + c2*Rkx*Rky - wB2*nx*ny;
  Lyz =           + c2*Rky*Rkz - wB2*ny*nz;
  Lzx =           + c2*Rkz*Rkx - wB2*nz*nx;

  % qubic equation for the frequncy: det(L)=0,
  % w^6 + a2*w^4 + a1*w^2 + a0 = 0
  a2 = Lxx + Lyy + Lzz - wL^2;
  a1 = Lxx*Lyy + Lyy*Lzz + Lzz*Lxx - Lxy^2 - Lyz^2 - Lzx^2 - wL^2*Lzz;
  a0 = Lxx*Lyy*Lzz + 2*Lxy*Lyz*Lzx - Lxx*Lyz^2 - Lyy*Lzx^2 - Lzz*Lxy^2;

  [n w1 w2 w3] = solve_cubic(1,a2,a1,a0);

  w1=sqrt(w1);
  w2=sqrt(w2);
  w3=sqrt(w3);
end

% function sutable for fortran
function [w1,w2,w3] = specfunc2(c1,c2,wL,wB2, kv, ak, bk, an, bn);

  % components of the n vector
  nx = sind(bn)*cosd(an);
  ny = sind(bn)*sind(an);
  nz = cosd(bn);

  % components of the k vector
  kx = kv*sind(bk)*cosd(ak);
  ky = kv*sind(bk)*sind(ak);
  kz = kv*cosd(bk);

  % rotated k vector, R_{aj} k_j
  ct=-1/4; st=sqrt(15)/4;
  kn  = kx*nx + ky*ny + kz*nz; % (k*n)
  kl = kz*ct + nz*kn*(1-ct) - (kx*ny-ky*nx)*st;

  % he3 parameters - not now
%  cpar2 = he3_cpar(ttc, P)^2;
%  cper2 = he3_cperp(ttc, P)^2;
%  wB2 = (2*pi*he3_nu_b(ttc, P))^2;
%  wL = he3_gyro*H;

  % qubic equation for the frequncy: det(L)=0,
  % w^6 + a2*w^4 + a1*w^2 + a0 = 0

  a2 = - (3*c1 - c2)*kv^2 - wB2 - wL^2;
  a1 = ...
    + kv^4*c1*(3*c1 - 2*c2) ...
    + kv^2*wB2*(2*c1 - c2) ...
    + kn^2*wB2*c2 ...
    + wL^2*(c1*kv^2 - c2*kl^2 + wB2*nz^2);
  a0 = ...
     + kv^4*c1*(kv^2*c1 + wB2)*(c2-c1) ...
     - wB2*kv^2*c1*c2*kn^2;

# 1. n_j k_j
# 2. |k|
# 3. n_z
# 4. l_j k_j

  [n w1 w2 w3] = solve_cubic(1,a2,a1,a0);

  w1=sqrt(w1);
  w2=sqrt(w2);
  w3=sqrt(w3);
end

% function sutable for fortran kx^2(w)
% only for bn=90deg!
function [k2a,k2b,k2c] = specfunc_k(c1,c2,wL,wB2, w, an, bn);

  % components of the n vector
  nx = sind(bn)*cosd(an);
  ny = sind(bn)*sind(an);
  nz = cosd(bn);

  % rotated k vector, R_{aj} k_j
  ct=-1/4; st=sqrt(15)/4;
  Rxx = ct + (1-ct)*nx*nx;
  Ryx = (1-ct)*ny*nx + st*nz;
  Rzx = (1-ct)*nz*nx - st*ny;

  Axx = - c1 + c2*Rxx*Rxx;
  Ayy = - c1 + c2*Ryx*Ryx;
  Azz = - c1 + c2*Rzx*Rzx;
  Axy =        c2*Rxx*Ryx;
  Ayz =        c2*Ryx*Rzx;
  Azx =        c2*Rzx*Rxx;
  Bxx = - wB2*nx*nx;
  Byy = - wB2*ny*ny;
  Bzz = - wB2*nz*nz;
  Bxy = - wB2*nx*ny;
  Byz = - wB2*ny*nz;
  Bzx = - wB2*nz*nx;

  % qubic equation for the frequncy: det(L)=0,
  % w^6 + a2*w^4 + a1*w^2 + a0 = 0
  a3 = Axx*Ayy*Azz + 2D0*Axy*Ayz*Azx ...
     - Axx*Ayz^2 - Ayy*Azx^2 - Azz*Axy^2;

  a2 = (Axx*Ayy + Ayy*Azz + Azz*Axx)*w^2 ...
     - (Axy^2 + Ayz^2 + Azx^2)*w^2 ...
     + Axx*Ayy*Bzz + Axx*Byy*Azz + Bxx*Ayy*Azz ...
     + 2D0*(Axy*Ayz*Bzx + Axy*Byz*Azx + Bxy*Ayz*Azx) ...
     - 2D0*(Axx*Ayz*Byz + Ayy*Azx*Bzx + Azz*Axy*Bxy) ...
     - Bxx*Ayz^2 - Byy*Azx^2 - Bzz*Axy^2;

   a1 = (Axx + Ayy + Azz)*w^4 - Azz*wL^2*w^2 ...
      + (Axx*Byy + Ayy*Bxx + Ayy*Bzz)*w^2 ...
      + (Azz*Byy + Azz*Bxx + Axx*Bzz)*w^2 ...
      - 2D0*(Axy*Bxy + Ayz*Byz + Azx*Bzx)*w^2 ...
      + Axx*Byy*Bzz + Bxx*Ayy*Bzz + Bxx*Byy*Azz ...
      + 2D0*(Axy*Byz*Bzx + Bxy*Ayz*Bzx + Bxy*Byz*Azx) ...
      - 2D0*(Bxx*Ayz*Byz + Byy*Azx*Bzx + Bzz*Axy*Bxy) ...
      - Axx*Byz^2 - Ayy*Bzx^2 - Azz*Bxy^2;

   a0 = w^6 + (Bxx+Byy+Bzz-wL^2)*w^4 ...
      + (Bxx*Byy + Byy*Bzz + Bzz*Bxx)*w^2 ...
      - (Bxy^2 + Byz^2 + Bzx^2 + Bzz*wL^2)*w^2 ...
      + Bxx*Byy*Bzz + 2D0*Bxy*Byz*Bzx ...
      - Bxx*Byz^2 - Byy*Bzx^2 - Bzz*Bxy^2;

  [n k2a k2b k2c] = solve_cubic(a3,a2,a1,a0);

end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> B phase magnon spectra.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! simple formula

      subroutine he3b_spec_simp(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        real*8 wL,wB2, cpar2t,cper2t,cpar2l,cper2l,c2l,c2t
        real*8 ct,st, kx,ky,kz, nx,ny,nz, kn,kl


        ! calculate (k*l)^2 = (h_a R_aj k_j)^2  for unit vectors
        ct = -0.25D0         ! \cos\theta
        st = sqrt(15D0/16D0) ! \sin\theta

        ! components of the n vector
        nx = dsin(bn)*dcos(an)
        ny = dsin(bn)*dsin(an)
        nz = dcos(bn)

        ! components of the unit k vector
        kx = dsin(bk)*dcos(ak)
        ky = dsin(bk)*dsin(ak)
        kz = dcos(bk)

        kn = kx*nx + ky*ny + kz*nz
        kl = kz*ct + kn*nz*(1D0-ct) - (kx*ny-ky*nx)*st


        ! Larmor and Leggett frequencies
        wL = he3_gyro*H
        wB2 = (const_2pi * he3_nu_b(ttc,P))**2


        ! c_par and c_perp for longitudinal and transverse magnons
        cpar2t = he3_cpar(ttc, P)**2
        cper2t = he3_cperp(ttc, P)**2
        cpar2l = he3_clpar(ttc, P)**2
        cper2l = he3_clperp(ttc, P)**2

        ! spin wave velocities for the given direction of k
        c2t = cper2t + (cpar2t-cper2t)*kl**2
        c2l = cper2l + (cpar2l-cper2l)*kl**2

        ! Acoustic magnon
        w1 = - wL/2D0 + sqrt((wL/2D0)**2
     .          + kv**2*c2t + 0.5D0*wB2*(1D0-nz**2))

        ! Optical magnon
        w2 = + wL/2D0 + sqrt((wL/2D0)**2
     .          + kv**2*c2t + 0.5D0*wB2*(1D0-nz**2))

        ! Longitudinal magnon
        w3 =  sqrt(kv**2*c2l + wB2*nz**2);

      end

!> He3-B acoustic magnon spectrum, simple formula
      function he3b_spec1s(ttc,P,H, kv, ak,bk, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_simp(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec1s = w1
      end

!> He3-B optical magnon spectrum, simple formula
      function he3b_spec2s(ttc,P,H, kv, ak,bk, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_simp(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec2s = w2
      end

!> He3-B longitudinal magnon spectrum, simple formula
      function he3b_spec3s(ttc,P,H, kv, ak,bk, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_simp(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec3s = w3
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! complete formula

      subroutine he3b_spec_full(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        real*8 wL,wB2, c1,c2
        real*8 ct,st, kx,ky,kz, nx,ny,nz, kn, Rkx,Rky,Rkz
        real*8 Lxx,Lyy,Lzz,Lxy,Lyz,Lzx,a0,a1,a2

        ! calculate (k*l)^2 = (h_a R_aj k_j)^2  for unit vectors
        ct = -0.25D0         ! \cos\theta
        st = sqrt(15D0/16D0) ! \sin\theta

        ! components of the n vector
        nx = dsin(bn)*dcos(an)
        ny = dsin(bn)*dsin(an)
        nz = dcos(bn)

        ! components of the unit k vector
        kx = dsin(bk)*dcos(ak)
        ky = dsin(bk)*dsin(ak)
        kz = dcos(bk)

        ! rotated unit k vector
        kn = kx*nx + ky*ny + kz*nz
        Rkx = kx*ct + nx*kn*(1-ct) - (ky*nz-kz*ny)*st
        Rky = ky*ct + ny*kn*(1-ct) - (kz*nx-kx*nz)*st
        Rkz = kz*ct + nz*kn*(1-ct) - (kx*ny-ky*nx)*st

        ! Larmor and Leggett frequencies
        wL = he3_gyro*H
        wB2 = (const_2pi * he3_nu_b(ttc,P))**2

        ! gradient terms
        c1 = kv**2 * he3_clperp(ttc,P)**2 ! K
        c2 = kv**2 * (he3_clperp(ttc,P)**2 - he3_clpar(ttc,P)**2) ! K'

        ! components of the L matrix (symmetric)
        Lxx = - c1 + c2*Rkx*Rkx - wB2*nx*nx
        Lyy = - c1 + c2*Rky*Rky - wB2*ny*ny
        Lzz = - c1 + c2*Rkz*Rkz - wB2*nz*nz
        Lxy =        c2*Rkx*Rky - wB2*nx*ny
        Lyz =        c2*Rky*Rkz - wB2*ny*nz
        Lzx =        c2*Rkz*Rkx - wB2*nz*nx

        ! qubic equation for the frequncy: det(L)=0,
        ! w^6 + a2*w^4 + a1*w^2 + a0 = 0
        a2 = Lxx + Lyy + Lzz - wL**2
        a1 = Lxx*Lyy + Lyy*Lzz + Lzz*Lxx
     .     - Lxy**2 - Lyz**2 - Lzx**2 - wL**2*Lzz
        a0 = Lxx*Lyy*Lzz + 2D0*Lxy*Lyz*Lzx
     .     - Lxx*Lyz**2 - Lyy*Lzx**2 - Lzz*Lxy**2

        call solve_cubic(1D0,a2,a1,a0, w1,w2,w3)
        w1=dsqrt(w1)
        w2=dsqrt(w2)
        w3=dsqrt(w3)
      end

!> He3-B magnon spectrum, full formula, lowest mode
      function he3b_spec1(ttc,P,H, kv, ak,bk, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_full(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec1 = w1
      end

!> He3-B magnon spectrum, full formula, middle mode
      function he3b_spec2(ttc,P,H, kv, ak,bk, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_full(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec2 = w2
      end

!> He3-B magnon spectrum, full formula, highest mode
      function he3b_spec3(ttc,P,H, kv, ak,bk, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_full(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec3 = w3
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inverse formula, k^2(w) for spin-wave propagating in x direction.
!
      subroutine he3b_spec_kx2(ttc,P,H, w, an,bn, k2a,k2b,k2c)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, w, an,bn, k2a,k2b,k2c
        real*8 ct,st, nx,ny,nz
        real*8 Rxx,Ryx,Rzx
        real*8 wL,wB2, c1,c2
        real*8 Axx,Ayy,Azz,Axy,Ayz,Azx
        real*8 Bxx,Byy,Bzz,Bxy,Byz,Bzx
       real*8 a0,a1,a2,a3

        ! theta is 104deg
        ct = -0.25D0         ! \cos\theta
        st = sqrt(15D0/16D0) ! \sin\theta

        ! components of the n vector
        nx = dsin(bn)*dcos(an)
        ny = dsin(bn)*dsin(an)
        nz = dcos(bn)

        ! components of the order parameter matrix
        Rxx = ct + (1-ct)*nx*nx
        Ryx = (1-ct)*ny*nx + st*nz
        Rzx = (1-ct)*nz*nx - st*ny

        ! Larmor and Leggett frequencies
        wL = he3_gyro*H
        wB2 = (const_2pi * he3_nu_b(ttc,P))**2

        ! gradient terms
        c1 = he3_clperp(ttc,P)**2 ! K
        c2 = (he3_clperp(ttc,P)**2 - he3_clpar(ttc,P)**2) ! K'

        ! components of the L = A*k^2 + B matrix
        Axx = - c1 + c2*Rxx*Rxx
        Ayy = - c1 + c2*Ryx*Ryx
        Azz = - c1 + c2*Rzx*Rzx
        Axy =        c2*Rxx*Ryx
        Ayz =        c2*Ryx*Rzx
        Azx =        c2*Rzx*Rxx
        Bxx = - wB2*nx*nx
        Byy = - wB2*ny*ny
        Bzz = - wB2*nz*nz
        Bxy = - wB2*nx*ny
        Byz = - wB2*ny*nz
        Bzx = - wB2*nz*nx

        ! Here I substituted L=Ak^2+B into equation for w
        ! and rewrite it as a qubic equation for k^2.

        a3 = Axx*Ayy*Azz + 2D0*Axy*Ayz*Azx
     .     - Axx*Ayz**2 - Ayy*Azx**2 - Azz*Axy**2

        a2 = (Axx*Ayy + Ayy*Azz + Azz*Axx)*w**2
     .     - (Axy**2 + Ayz**2 + Azx**2)*w**2
     .    + Axx*Ayy*Bzz + Axx*Byy*Azz + Bxx*Ayy*Azz
     .    + 2D0*(Axy*Ayz*Bzx + Axy*Byz*Azx + Bxy*Ayz*Azx)
     .    - 2D0*(Axx*Ayz*Byz + Ayy*Azx*Bzx + Azz*Axy*Bxy)
     .    - Bxx*Ayz**2 - Byy*Azx**2 - Bzz*Axy**2

        a1 = (Axx + Ayy + Azz)*w**4 - Azz*wL**2*w**2
     .     + (Axx*Byy + Ayy*Bxx + Ayy*Bzz)*w**2
     .     + (Azz*Byy + Azz*Bxx + Axx*Bzz)*w**2
     .     - 2D0*(Axy*Bxy + Ayz*Byz + Azx*Bzx)*w**2
     .     + Axx*Byy*Bzz + Bxx*Ayy*Bzz + Bxx*Byy*Azz
     .     + 2D0*(Axy*Byz*Bzx + Bxy*Ayz*Bzx + Bxy*Byz*Azx)
     .     - 2D0*(Bxx*Ayz*Byz + Byy*Azx*Bzx + Bzz*Axy*Bxy)
     .     - Axx*Byz**2 - Ayy*Bzx**2 - Azz*Bxy**2

        a0 = w**6 + (Bxx+Byy+Bzz-wL**2)*w**4
     .     + (Bxx*Byy + Byy*Bzz + Bzz*Bxx)*w**2
     .     - (Bxy**2 + Byz**2 + Bzx**2 + Bzz*wL**2)*w**2
     .     + Bxx*Byy*Bzz + 2D0*Bxy*Byz*Bzx
     .     - Bxx*Byz**2 - Byy*Bzx**2 - Bzz*Bxy**2

        call solve_cubic(a3,a2,a1,a0, k2a,k2b,k2c)
      end

!> He3-B magnon spectrum, kx^2(w), lowest mode
      function he3b_spec_kx2a(ttc,P,H, w, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, w, an,bn, k2a,k2b,k2c
        call he3b_spec_kx2(ttc,P,H, w, an,bn, k2a,k2b,k2c)
        he3b_spec_kx2a = k2a
      end

!> He3-B magnon spectrum, kx^2(w), middle mode
      function he3b_spec_kx2b(ttc,P,H, w, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, w, an,bn, k2a,k2b,k2c
        call he3b_spec_kx2(ttc,P,H, w, an,bn, k2a,k2b,k2c)
        he3b_spec_kx2b = k2b
      end

!> He3-B magnon spectrum, kx^2(w), highest mode
      function he3b_spec_kx2c(ttc,P,H, w, an,bn) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, w, an,bn, k2a,k2b,k2c
        call he3b_spec_kx2(ttc,P,H, w, an,bn, k2a,k2b,k2c)
        he3b_spec_kx2c = k2c
      end

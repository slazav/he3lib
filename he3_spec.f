!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! B phase magnon spectra.

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

      function he3b_spec1s(ttc,P,H, kv, ak,bk, an,bn)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_simp(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec1s = w1
      end

      function he3b_spec2s(ttc,P,H, kv, ak,bk, an,bn)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_simp(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec2s = w2
      end

      function he3b_spec3s(ttc,P,H, kv, ak,bk, an,bn)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_simp(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec3s = w3
      end



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

      function he3b_spec1(ttc,P,H, kv, ak,bk, an,bn)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_full(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec1 = w1
      end

      function he3b_spec2(ttc,P,H, kv, ak,bk, an,bn)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_full(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec2 = w2
      end

      function he3b_spec3(ttc,P,H, kv, ak,bk, an,bn)
        implicit none
        include 'he3.fh'
        real*8 ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3
        call he3b_spec_full(ttc,P,H, kv, ak,bk,an,bn, w1,w2,w3)
        he3b_spec3 = w3
      end


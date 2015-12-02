! B phase in a strong magnetic field (B2 phase).

! Calculation of B-phase gap distortion and spin polarization.
! Based on Ashida and Nagai paper (Progr.Theor.Phys. 74 949 (1985)).

! Free energy derivatives: F1=dF/dGap1^2,  F2=dF/dGap1^2, F3=dF/we
      subroutine he3_b2_fder(ttc,gap, A,B,we,w0,f0a, F1,F2,F3)
        implicit none
        real*8 ttc ! T/Tc
        real*8 gap ! Gap0/Tc
        real*8 A,B ! gap distortions, we/w0
        real*8 we,w0    ! effective and external fields in ttc units
        real*8 f0a      ! fermi-liquid parameter
        real*8 F1,F2,F3 ! returned values

        ! integration is done exactly as in he3_math.f:math_dint2d function,
        ! but we want to calculate three integrals at once.
        ! Integration region: [0:1]x[0:1] square (1/(xi-1) and \cos\theta).
        ! Integrand is good and we do not need any smart integration with lots of points

        real*8 xi,dxi, g1,g2,g0,gg,gz
        real*8 Ep,Em,Ez,Eb, Ap,Am,Ab, I1,I2,I3
        real*8 ss(2), dx, dy, x,y
        integer ix, iy, nx, ny, ixx,iyy
        nx=30
        ny=30
        ss(1) = 0.5D0 - 0.5D0/dsqrt(3D0)
        ss(2) = 0.5D0 + 0.5D0/dsqrt(3D0)
        dx=1D0/dble(nx)
        dy=1D0/dble(ny)
        F1=0D0
        F2=0D0
        F3=0D0

        g0 = gap**2
        g1 = (gap*(1D0+A))**2
        g2 = (gap*(1D0+B))**2

        do ix=1,nx    ! x coordinate is 1/(xi-1)
          do iy=1,ny  ! y coordinate is pz = \cos\theta
            do ixx=1,2
              do iyy=1,2
                x = dx*(dble(ix) - ss(ixx))
                y = dy*(dble(iy) - ss(iyy))
                xi = 1D0/x-1D0
                dxi = -1D0/x**2  ! dxi/dx

                gg = g1*(1D0-y**2) + g2*y**2
                Ez = dsqrt(xi**2 + g2*y**2);
                Eb = dsqrt(xi**2 + g0);
                Ep = dsqrt(xi**2 + gg + we**2/4D0 - we*Ez)
                Em = dsqrt(xi**2 + gg + we**2/4D0 + we*Ez)

                ! A = (f(E)-1/2)/2E = dF/dE*dE/d(sqrt(E))
                Ab = -0.25D0*dtanh(Eb/(2D0*ttc))/Eb
                Ap = -0.25D0*dtanh(Ep/(2D0*ttc))/Ep
                Am = -0.25D0*dtanh(Em/(2D0*ttc))/Em

                ! integrands
                F1 = F1 + (Ap+Am-2D0*Ab)*(1D0-y**2) *dxi*dx*dy/4D0
                F2 = F2 + (Ap*(1D0-we/(2D0*Ez))
     .                  +  Am*(1D0+we/(2D0*Ez))
     .                  -  2D0*AB)*y**2 *dxi*dx*dy/4D0
                F3 = F3 + (Ap*(we/2D0-Ez)
     .                  +  Am*(we/2D0+Ez)) *dxi*dx*dy/4D0
              enddo
            enddo
          enddo
        enddo
        F3 = F3 + (we-w0)/(4D0*f0a)
      end

! Minimum of the free energy F(gap1,gap2,He)
! Returns gap1,gap2,He
      subroutine he3_b2_fmin(ttc,p,H, gap1,gap2,He)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,H ! T/Tc; pressure, bar; field, G
        real*8 gap1,gap2,He
        real*8 A,B,we  ! gap distortions, effective field

        real*8 w0,gap,f0a,tc
        real*8 F1a,F2a,F3a, F1b,F2b,F3b, sl1,sl2,sl3, dA,dB,dwe
        integer i

        gap = he3_gap(ttc,p)
        f0a = he3_f0a(p)
        ! f0a = -0.75D0 ! for tests: this value was used in the original paper
        tc  = he3_tc(p)*1D-3*const_kb   ! Tc in erg units
        w0  = H*he3_gyro*const_hbar/tc  ! w0 in tc units

        A=0D0
        B=0D0
        we=2D0*w0
        call he3_b2_fder(ttc,gap,A,B,we,w0,f0a, F1a,F2a,F3a)

        sl1 = 3D0 ! slope
        sl2 = 3D0 ! slope
        sl3 = 3D0 ! slope
        do i=1,100
          dA  = sl1*F1a
          dB  = sl2*F2a
          dwe = sl3*F3a
          A   = A + dA
          B   = B + dB
          we  = we + dwe
          call he3_b2_fder(ttc,gap,A,B,we,w0,f0a, F1b,F2b,F3b)
          if ((abs(F1b)+abs(F2b)+abs(F3b)).lt.1D-10) goto 1
          if (A.le.-1D0.or.B.le.-1D0) then
            A=1D0/0D0
            B=1D0/0D0
            we=1D0/0D0
            goto 1
          endif
          if (abs(dA/(F1a-F1b))<1D2)  sl1 = (dA/(F1a-F1b))
          if (abs(dB/(F2a-F2b))<1D2)  sl2 = (dB/(F2a-F2b))
          if (abs(dwe/(F3a-F3b))<1D2) sl3 = (dwe/(F3a-F3b))
          F1a=F1b
          F2a=F2b
          F3a=F3b
        enddo
1       gap1=gap*(1D0+A)
        gap2=gap*(1D0+B)
        He = we*tc/(he3_gyro*const_hbar) ! Tc units -> G
      end

      function he3_b2gap1(ttc,p,H)
        implicit none
        real*8 ttc,p,H
        real*8 gap1,gap2,He, he3_b2gap1
        call he3_b2_fmin(ttc,p,H, gap1,gap2,He)
        he3_b2gap1 = gap1
      end

      function he3_b2gap2(ttc,p,H)
        implicit none
        real*8 ttc,p,H
        real*8 gap1,gap2,He, he3_b2gap2
        call he3_b2_fmin(ttc,p,H, gap1,gap2,He)
        he3_b2gap2 = gap2
      end

      function he3_b2heff(ttc,p,H)
        implicit none
        real*8 ttc,p,H
        real*8 gap1,gap2,He, he3_b2heff
        call he3_b2_fmin(ttc,p,H, gap1,gap2,He)
        he3_b2heff = He
      end

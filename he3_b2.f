! B phase in a strong magnetic field (B2 phase).

!! B_ab [mk] vs P, ttc
! Inseob Hahn PhD thesis, p79
! see also code at http://spindry.phys.northwestern.edu/he3.htm
      function pfunc(p,a,n)
      ! function used to approximate pressure-dependent parameters
        real*8 p,a(n),pfunc
        integer i
        pfunc=0
        do i=1,n-1
          pfunc = pfunc + a(i)*p**(i-1)
        enddo
        pfunc = pfunc/(1D0+a(n)*p)
      end
      function He3_b2hcr(ttc,P)
        implicit none
        include 'he3.fh'
        real*8 P,ttc, tc,tab
        real*8 p0,p1,p2,q1,pp
        real*8 f1,f2,f3,f4,f5
        real*8 Bc,Bo,gr, Pa, BB, X2
        real*8 f0a,pfunc

        ! values for pfunc()
        real*8 coeff_Bc(4), coeff_Ta(4), coeff_f5(6)
        data
     .    coeff_Bc /3391D0, 21500D0, -8490D0, 2.098D0/,
     .    coeff_Ta /1.5745D0,-1.1222D0,0.3242D0,0D0/,
     .    coeff_f5 /0.041870D0, 5.417531D0, -10.044312D0,
     .              7.639438D0, -2.379517D0, -0.537422D0/

        Pa = 34.338
        pp = P/Pa

        f0a = he3_f0a(P)
        Tc  = he3_tc(P)

       ! Critical field at zero temperature.
        Bc  = pfunc(pp, coeff_Bc, 4)

       ! Tab/Tc. We need this value which is defined
       ! in the whole pressure range and goes above 1 at low pressures.
       ! This value is close to he3_tab(P)/he3_tc(P)
        tab = pfunc(pp, coeff_Ta, 4)

        f3 = 1.41D0
        f4 = -0.29D0 -0.41D0*pp
        gr = 0.616D0 - 1.174D0*pp + 0.301D0*pp**2

        Bo = 19700 * Tc * (1+f0a)
        BB = (Bo/Bc)**2 *gr/4D0
        ! This is the original value of f5 which has problems near Tab=1
        ! In expression A/BB both parts crosses 0 at this point.
        ! And this crossing comes from different places (gr and tab definition)
        !f5 = (1D0 - f3*Tab**6 - f4*Tab**8 - (1D0-f3-f4)*Tab**2
     .  !      + (1D0+2D0*f3+3D0*f4)*(Tab**4-Tab**2))
     .  !     /(BB * (Tab**4-Tab**2)) - 1D0
        ! I will use the smooth approximation instead:
        f5 = pfunc(pp, coeff_f5, 6)

        f2 = BB * (1D0+f5) - (1D0+2D0*f3+3D0*f4)
        f1 = 1D0+f5-f2-f3-f4

        X2 = (f1*ttc**2 + f2*ttc**4 + f3*ttc**6 + f4*ttc**8)
     .    / (1+f5*ttc**2)
        he3_b2hcr = dsqrt(1D0 - X2**2) * Bc
      end

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

        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          gap1=1D0/0D0
          gap2=1D0/0D0
          He=1D0/0D0
          return
        endif
        if (ttc.eq.1D0) then
          gap1=0
          gap2=0
          He=H/(1D0+he3_f0a(p))
          return
        endif

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
        do i=1,1000
          dA  = sl1*F1a
          dB  = sl2*F2a
          dwe = sl3*F3a
          A   = A + dA
          B   = B + dB
          we  = we + dwe
          call he3_b2_fder(ttc,gap,A,B,we,w0,f0a, F1b,F2b,F3b)
          if ((abs(F1b)+abs(F2b)+abs(F3b)).lt.1D-10) goto 2
          if (A.le.-1D0.or.B.le.-1D0) goto 1
          if (abs(dA/(F1a-F1b))<1D2)  sl1 = (sl1 + dA/(F1a-F1b))/2D0
          if (abs(dB/(F2a-F2b))<1D2)  sl2 = (sl1 + dB/(F2a-F2b))/2D0
          if (abs(dwe/(F3a-F3b))<1D2) sl3 = (sl1 + dwe/(F3a-F3b))/2D0
          F1a=F1b
          F2a=F2b
          F3a=F3b
        enddo
1       A=1D0/0D0
        B=1D0/0D0
        we=1D0/0D0
2       gap1=gap*(1D0+A)
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

!HH> Vibrating wire calibration programs from Lancaster ULT (original version)

!> These functions produce same results as old Lancaster code, I do not plan to change them.
!> I plan to do a separate version using functions from the library (Yosida functions, viscosity, etc.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Calibration of vibrating wire in mixing chamber (diluted phase), modified, complex: freq + 1i*width
!> Arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_newwire_d(t, rho, diam, fre) !FC>
        implicit none
        include 'he3.fh'
        real*8 t, rho, diam, fre
        real*8 alpha, ee, rad, vol, vol3, mstar
        real*8 conc, cratio, rho3, rratio, eta, pend
        real*8 G,L, b, k,kp, k2,k3
        complex*16 kk
        real*8 ff,df

        alpha = 2.2D0  ! mfp fudge
        ee = 10D0      ! switch strength to ballistic
        rad = diam/2e6 ! radius in m
        vol = 27.58D0  ! he4 molar volume
        mstar = 2.46D0 !

        ! he3 concentration
        conc = 0.066D0 + 0.5056D0*t**2
     .       - 0.2488D0*t**3 + 18.22D0*t**4
     .       - 74.22D0*t**5
        cratio=0.0665/conc

        vol3 = vol*(1D0 + 0.286D0*conc)
        rho3 = conc * 3.016D0*mstar/vol3
        rratio = rho3/rho;

        ! viscosity
        if (t.lt.0.0165D0) then
          eta=0.305D-7/t**2 + 1.35D-7/t + 2.2D-6
        elseif (t.gt.0.068D0) then
          eta=0.29D-7/t**2 + 1.65D-7/t + 2.3D-6
          ! improved version, allowing for mfp effects, reasonable fit up to 100mK
        else
          eta = 0.277D-7/t**2 + 3.4D-7/t  !original Zeegers et al
        endif

        ! viscous penetration depth
        pend = dsqrt(eta/(2000D0*const_pi*rho3*fre))

        ! Stokes G and mfp L
        G = rad/pend;
        L = eta/105.5302 * cratio**(4D0/3D0)

        kk = math_stokes(G)
        k  = real(kk)
        kp = imag(kk)
        b = 0.25D0 * 0.579D0 * L/rad
        b = b*(1D0 + ee*alpha*L/rad) / (1D0 + ee*L/rad)

        k2 = 1D0 + (k - 1D0) /
     .       ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

        k3 = (kp + G**2*b*((k-1D0)**2 + kp**2)) /
     .       ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

        df = fre * rratio * k3 * (1D0-1.14D0*rratio*k2)
        ff = fre * rratio * 0.5D0 * k2 * (1D0-0.75*rratio*k2)
        he3_newwire_d = dcmplx(ff, df)
      end

!> Calibration of vibrating wire in mixing chamber (diluted phase), frequency
!> arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_newwire_d_f(t, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, rho, diam, fre
        he3_newwire_d_f = real(he3_newwire_d(t, rho, diam, fre))
      end

!> Calibration of vibrating wire in mixing chamber (diluted phase), width
!> arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_newwire_d_w(t, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, rho, diam, fre
        he3_newwire_d_w = imag(he3_newwire_d(t, rho, diam, fre))
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Calibration of vibrating wire in superfluid He3-B (modified), complex: freq + 1i*width
!> arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz].
      function he3_newwire_b(t, p, rho, diam, fre) !FC>
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        real*8 alpha, ee, rad
        real*8 tc, vol, rho3, rratio, f1s
        real*8 vistc,vis, eovl
        real*8 ttc, gap, y0,y1,y2,y3,y5,y6, ts
        real*8 pend, zeta
        real*8 G,L, k,kp, b, k2,k3
        complex*16 kk
        real*8 ff,df
        real*8 lanc_redvis, lanc_visc

        alpha=1.9D0     ! mfp fudge
        rad = diam/2D6  ! radius in m
        ee  = 1D0       ! ballistic switch wick

        tc    = he3_tc(p)/1000 ! Tc in K
        vol   = he3_vm(p)      ! Molar volume
        f1s = he3_f1s(p);

        eovl = 0.2D0*(6.023D29/vol)**(4D0/3D0)
     .       * (3D0*9.8696D0)**(1D0/3D0)*1.0546D-34

        vistc=lanc_visc(tc,p)
        ttc = t/tc    !reduced temperature
!        L = he3_visc_fpath(ttc,p) !mean free path

        vis=vistc*lanc_redvis(ttc) ! fudged effective viscosity from CHH

        gap = he3_gap(ttc,p)
        y0 = he3_yosida(ttc, gap, 0)
        y1 = he3_yosida(ttc, gap, 1)
        y2 = he3_yosida(ttc, gap, 2)
        y3 = he3_yosida(ttc, gap, 3)

        y5 = (8D0/15D0*y2/y1 + 5D0/8D0*y3/y2)/y2
        y6 = dsqrt(y2*y0)

        rho3  = he3_rho_nb(ttc,p)*he3_rho(p) ! effective density
        rratio = rho3/rho                    ! density ratio

        pend = dsqrt(vis/(2000D0*const_pi*rho3*fre)) ! penetration depth
        zeta=0.5D0*y5*vis/eovl                       ! effective slip length
        L=vis/(eovl*y6)                         ! mean free path
        alpha = 1.156D0*alpha/(y5*y6)                ! effective alpha
        G = rad/pend                                 ! gamma for Stokes

        kk=math_stokes(G)
        k  = real(kk)
        kp = imag(kk)
        b = 0.25D0 * zeta/rad
        b = b*(1D0 + ee*alpha*L/rad) / (1D0 + ee*L/rad)

        k2 = 1D0 + (k - 1D0) /
     .       ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

        k3 = (kp + G**2*b*((k-1D0)**2 + kp**2)) /
     .       ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

        df = fre * rratio * k3 ! * (1D0-1.14D0*rratio*k2)
        ff = fre * rratio * 0.5D0 * k2 ! * (1D0-0.75*rratio*k2)

        he3_newwire_b = dcmplx(ff, df)
      end

!> Calibration of vibrating wire in superfluid He3-B, frequency
!> arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_newwire_b_f(t, p, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        he3_newwire_b_f =
     .    real(he3_newwire_b(t, p, rho, diam, fre))
      end

!> Calibration of vibrating wire in superfluid He3-B, width
!> arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_newwire_b_w(t, p, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        he3_newwire_b_w =
     .    imag(he3_newwire_b(t, p, rho, diam, fre))
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Calibration of vibrating wire in normal He3 (modified), complex: freq + 1i*width
!> Cylinder programme using wide line treatment and slip fudge.
!> arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz].
      function he3_newwire_n(t, p, rho, diam, fre) !FC>
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        real*8 alpha,rad,ee,vol,eovl,vis,rho_h,rhorat,pen
        real*8 G,L,b,k,kp,k2,k3
        complex*16 kk
        real*8 ff,df,sq1,sq2
        real*8 lanc_visc
        integer j

        alpha=1.9D0     ! mfp fudge
        rad = diam/2D6  ! radius in m
        ee  = 1D0       ! ballistic switch

        vol=he3_vm(p)   ! molar volume

        eovl=0.2D0*(6.023D29/vol)**(4D0/3D0)*
     .             (3D0*9.8696D0)**(1D0/3D0)*1.0546D-34

        vis=lanc_visc(t,p)
        rho_h=3.016D0/vol
        rhorat=rho_h/rho

        !! shift using approx formula for fork; remove this section for wires
        !  x=dlog10(t*1000);
        !  F1a=1990D0 - 2165D0*x + 1150D0*x**2 - 247D0*x**3 + 8.65D0*x**4 +2.57D0*x**5
        !  f=fre-F1a

        ff=fre

        do j=1,6
          pen=dsqrt(vis/(2000D0*const_pi*rho_h*ff))
          G=rad/pen
          L=vis/eovl
          kk=math_stokes(G)
          k  = real(kk)
          kp = imag(kk)

          b = 0.25D0 * 0.579D0 * L/rad
          b = b*(1D0 + ee*alpha*L/rad) / (1D0 + ee*L/rad)

          k2 = 1D0 + (k - 1D0) /
     .         ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

          k3 = (kp + G**2*b*((k-1D0)**2 + kp**2)) /
     .         ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

          ff=fre*(dsqrt(1D0+(rhorat*k2/2D0)**2)-rhorat*k2/2D0)
        enddo

        df=rhorat*k3*ff
        sq1=ff*dsqrt((rhorat*(k3+k2))**2 + 4D0)
        sq2=ff*dsqrt((rhorat*(k3-k2))**2 + 4D0)

        df= df + sq1 - sq2  !exact width
        he3_newwire_n = dcmplx(ff, df)

      end

!> Calibration of vibrating wire in normal He3 (modified), freq
!> Cylinder programme using wide line treatment and slip fudge.
!> arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_newwire_n_f(t, p, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        he3_newwire_n_f =
     .    real(he3_newwire_n(t, p, rho, diam, fre))
      end

!> Calibration of vibrating wire in normal He3 (modified), width
!> Cylinder programme using wide line treatment and slip fudge.
!> arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_newwire_n_w(t, p, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        he3_newwire_n_w =
     .    imag(he3_newwire_n(t, p, rho, diam, fre))
      end




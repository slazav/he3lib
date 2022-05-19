!H> Vibrating wire calibration programs from Lancaster ULT (original version)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Mixing chamber calibration, mixture (original Lancaster version)
!arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
!output:    complex version f + i*df
      function he3_wire_mix_c(t, rho, diam, fre)
        implicit none
        include 'he3.fh'
        real*8 t, rho, diam, fre
        real*8 alpha, ee, rad, vol, vol3, mstar
        real*8 conc, cratio, rho3, rratio, eta, pend
        complex*16 he3_wire_mix_c
        real*8 G,L, k,kp, b, k2,k3
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

        call math_stokes(G,k,kp)
        b = 0.25D0 * 0.579D0 * L/rad
        b = b*(1D0 + ee*alpha*L/rad) / (1D0 + ee*L/rad)

        k2 = 1D0 + (k - 1D0) /
     .       ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

        k3 = (kp + G**2*b*((k-1D0)**2 + kp**2)) /
     .       ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

        df = fre * rratio * k3 * (1D0-1.14D0*rratio*k2)
        ff = fre * rratio * 0.5D0 * k2 * (1D0-0.75*rratio*k2)
        he3_wire_mix_c = dcmplx(ff, df)
      end

!> Mixing chamber calibration, frequency vs temperatere
!> arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_wire_mix_f(t, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, rho, diam, fre
        complex*16 he3_wire_mix_c
        he3_wire_mix_f = real(he3_wire_mix_c(t, rho, diam, fre))
      end

!> Mixing chamber calibration, width vs temperatere
!> arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_wire_mix_w(t, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, rho, diam, fre
        complex*16 he3_wire_mix_c
        he3_wire_mix_w = imag(he3_wire_mix_c(t, rho, diam, fre))
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! a plot of CHH  reduced viscosity gives the following reasonable form.
! Basically, LOG(redvis) vz sqrt(1-T)^.5 reasonably slow
      function redvis(ttc)
        implicit none
        real*8 ttc, redvis, SLOPE, D0, N0

        if (ttc.lt.0.6D0) then
          redvis=0.11D0
        elseif (ttc.lt.0.7D0) then
          SLOPE=-0.8562D0
          D0=0.4D0
          N0=-0.9586D0
          redvis=10D0**(N0 + SLOPE*(sqrt(1-ttc)-sqrt(D0)))
        elseif (ttc.lt.0.8D0) then
          SLOPE=-0.6183D0
          D0=0.3D0
          N0=-0.8861D0
          redvis=10D0**(N0+SLOPE*(sqrt(1-ttc)-sqrt(D0)))
        elseif (ttc.lt.0.9D0) then
          SLOPE=-1.4172D0
          D0=0.2D0
          N0=-0.8239D0
          redvis=10D0**(N0+SLOPE*(sqrt(1-ttc)-sqrt(D0)))
        elseif (ttc.lt.0.95D0) then
          SLOPE=-1.7352D0
          D0=0.1D0
          N0=-0.6383D0
          redvis=10D0**(N0+SLOPE*(sqrt(1-ttc)-sqrt(D0)))
        elseif (ttc.lt.0.975D0) then
          SLOPE=-1.6177D0
          D0=0.05D0
          N0=-0.4776D0
          redvis=10D0**(N0+SLOPE*(sqrt(1-ttc)-sqrt(D0)))
        elseif (ttc.lt.1D0) then
          SLOPE=-2.3503D0
          D0=0.025D0
          N0=-0.3716D0
          redvis=10D0**(N0+SLOPE*(sqrt(1-ttc)-sqrt(D0)))
        else
          redvis = 1D0
        endif
      end

      function visc(t, p)
        implicit none
        real*8 t, p, visc
        real*8 visca, viscb

        ! normal state viscosity
        !     eta= 1/( AT^2 + B)\
        !  eta in Pa s T in mK originally, now K
        !  using the data given in Carless, Hall and Hook JLTP 50,583,83\
        !  table 1 , page 593 as smoothed by AMG.\
        !  The data, converted to a Greywall T scale\
        !  by multiplication of T by 1.12.\
        if (p.lt.1.28D0) then
          visca=0.38D0 - 0.007D0*(1.28D0-p)/1.18D0
        elseif (p.lt.4.65D0) then
          visca=0.424D0 - 0.044D0*(4.65D0-p)/3.37D0
        elseif (p.lt.9.89) then
          visca=0.495D0 - 0.071D0*(9.89D0-p)/5.24D0
        elseif (p.lt.19.89) then
          visca=0.603D0 - 0.108D0*(19.89D0-p)/10D0
        else
          visca=0.603D0 + 0.107D0*(p-19.89D0)/9.45D0
        endif
        visca=visca * 10D0*1.12D0*1.12D0
        visca=visca*1D6  ! convert back to T in K

      ! normal state viscosity
      !
      ! eta= 1/( AT^2 + B)\
      !
      ! using the data given in Carless, Hall and Hook JLTP 50,583,83\
      ! table 1 , page 593 as smoothed by AMG.\
      ! The data, converted to a Greywall T scale\
      ! by multiplication of T by 1.12.\

        if (p.lt.1.28D0) then
          viscb = 0.06D0 - 0.02D0*(1.28D0-p)/1.18D0
        elseif (p.lt.4.65D0) then
          viscb = 0.19D0 - 0.13D0*(4.65D0-p)/3.37D0
        elseif (p.lt.9.89D0) then
          viscb = 0.43D0 - 0.24D0*(9.89D0-p)/5.24D0
        elseif (p.lt.19.89D0) then
          viscb = 0.94D0 - 0.56D0*(19.89D0-p)/10D0
        else
          viscb = 0.94D0 + 0.56D0*(p-19.89D0)/9.45D0
        endif
        viscb = viscb*10D0

        visc=1D0/(visca*t**2 + viscb)
      end

! Superfluid He3-B calibration (original Lancaster version)
! arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
! complex version f + i*df
      function he3_wire_bphase_c(t, p, rho, diam, fre)
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        real*8 alpha, ee, rad
        real*8 tc, vol, rho3, rratio, f1s
        real*8 visc,vistc,vis, eovl
        real*8 ttc, gap, y0,y1,y2,y3,y5,y6, ts
        real*8 pend, zeta
        complex*16 he3_wire_bphase_c
        real*8 G,L, k,kp, b, k2,k3
        real*8 ff,df
        real*8 redvis

        alpha=1.9D0     ! mfp fudge
        rad = diam/2D6  ! radius in m
        ee  = 1D0       ! ballistic switch wick

        tc    = he3_tc(p)/1000 ! Tc in K
        vol   = he3_vm(p)      ! Molar volume
        f1s = he3_f1s(p);

        eovl = 0.2D0*(6.023D29/vol)**(4D0/3D0)
     .       * (3D0*9.8696D0)**(1D0/3D0)*1.0546D-34

        vistc=visc(tc,p)
        ttc = t/tc    !reduced temperature
!        L = he3_visc_fpath(ttc,p) !mean free path

        vis=vistc*redvis(ttc) ! fudged effective viscosity from CHH

!        gap = he3_gap(ttc,p)
!        y0 = he3_yosida(ttc, gap, 0)
!        y1 = he3_yosida(ttc, gap, 1)
!        y2 = he3_yosida(ttc, gap, 2)
!        y3 = he3_yosida(ttc, gap, 3)
!
!        y5 = (8D0/15D0*y2/y1 + 5D0/8D0*y3/y2)/y2
!        y6 = dsqrt(y2*y0)


        ! Yosida functions. These formulase are very inaccurate,
        ! but we want to reproduce original calculation here
        ts = ttc*(0.9074D0 - 0.0075D0*ttc
     .        - 0.0216D0*ttc**2 + 0.1396D0*ttc**3 - 0.0611*ttc**4)

        if (ttc.lt.0.94D0) then
          y0 = dexp(-1.76388D0/ts)
     .     *(3.454D0 - 0.88D0*ts
     .       + 4.625D0*ts**2 - 1.367D0*ts**3)/dsqrt(ts)
        else
          y0=1.985D0*ts -0.985D0
        endif

        if (ttc.lt.0.80D0) then
          y5 = dexp(1.76388D0/ts)*(0.10177D0 + 1.1958D0*ts
     .      - 1.425D0*ts**2 + 0.392*ts**3)/dsqrt(ts)
        else
          y5 = dexp(1.76388D0/ts)*(0.19847D0
     .       + 0.335D0*sqrt(1-ts))/dsqrt(ts)
        endif

        if (ttc.lt.0.90D0) then
          y6 = exp(-1.76388D0/ts)*(2.402D0 + 0.4467D0*ts
     .        - 2.117D0*ts**2 + 4.1D0*ts**3)
        else
          y6=1D0 - sqrt(1D0-ts)*(-4.517D0
     .       + 13.275D0*ts - 7.5D0*ts**2)
        endif

        rho3 = (1D0+f1s/3D0)*y0/(1+f1s*y0/3D0)*he3_rho(p)
!        rho3  = he3_rho_nb(ttc,p)*he3_rho(p) ! effective density
        rratio = rho3/rho         ! density ratio

        pend = dsqrt(vis/(2000D0*const_pi*rho3*fre)) ! penetration depth
        zeta=0.5D0*y5*vis/eovl                       ! effective slip length
        L=vis/(eovl*y6)                         ! mean free path
        alpha = 1.156D0*alpha/(y5*y6)                ! effective alpha
        G = rad/pend                                 ! gamma for Stokes

        call math_stokes(G,k,kp)
        b = 0.25D0 * zeta/rad
        b = b*(1D0 + ee*alpha*L/rad) / (1D0 + ee*L/rad)

        k2 = 1D0 + (k - 1D0) /
     .       ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

        k3 = (kp + G**2*b*((k-1D0)**2 + kp**2)) /
     .       ((1D0 + G**2*b*kp)**2 + G**4*b**2*(k-1D0)**2)

        df = fre * rratio * k3 ! * (1D0-1.14D0*rratio*k2)
        ff = fre * rratio * 0.5D0 * k2 ! * (1D0-0.75*rratio*k2)

        he3_wire_bphase_c = dcmplx(ff, df)
      end

!> Superfluid He3-B calibration, frequency vs temperatere
!> arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_wire_bphase_f(t, p, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        complex*16 he3_wire_bphase_c
        he3_wire_bphase_f =
     .    real(he3_wire_bphase_c(t, p, rho, diam, fre))
      end

!> Superfluid He3-B calibration, width vs temperatere
!> arguments: temperature [K], pressure [bar], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_wire_bphase_w(t, p, rho, diam, fre) !F>
        implicit none
        include 'he3.fh'
        real*8 t, p, rho, diam, fre
        complex*16 he3_wire_bphase_c
        he3_wire_bphase_w =
     .    imag(he3_wire_bphase_c(t, p, rho, diam, fre))
      end


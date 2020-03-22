! Original wire calibration programs from Lancaster ULT

!     Mixing chamber calibration, mixture
!     arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
!     complex version f + i*df
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

!     Mixing chamber calibration, frequency vs temperatere
!     arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_wire_mix_f(t, rho, diam, fre)
        implicit none
        include 'he3.fh'
        real*8 t, rho, diam, fre
        complex*16 he3_wire_mix_c
        he3_wire_mix_f = real(he3_wire_mix_c(t, rho, diam, fre))
      end

!     Mixing chamber calibration, width vs temperatere
!     arguments: temperature [K], rho wire [g/cm^3], wire diameter [um], frequency [Hz]
      function he3_wire_mix_w(t, rho, diam, fre)
        implicit none
        include 'he3.fh'
        real*8 t, rho, diam, fre
        complex*16 he3_wire_mix_c
        he3_wire_mix_w = imag(he3_wire_mix_c(t, rho, diam, fre))
      end

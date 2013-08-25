! B phase spin diffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Collision integral for Bogoliubov quasiparticles

! Collision integral in Einzel approximation
! Einzel, Wolfle, Hirschfeld, JLTP80 (1990), Appendix, p.66
! + my small fixes
      function he3_coll_int(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 xi, ttc, gap, x
        real*8 J0,J1,J2,J3, K0,K1,K2,K3, I0,I1,I2,I3
        real*8 a0,a1,a2,a3, b0,b1,b2, g0,d0

        x = gap/ttc
        a0 = -0.5768578D0
        a1 =  0.2694D0
        a2 =  0.29D0
        a3 = -0.08D0
        J0 = 3D0/4D0 / const_pi**0.5D0
     .       * (1D0+2D0*x)**1.5D0
     .       / (dexp(x) + a0 + a1*x + a2*x**2 + a3*x**3)
        J1 = x**2 / (2D0*const_pi)**0.5D0
     .       * (0.5D0 + x)**1.5D0
     .       / (1.3D0 + x**2) * dexp(-x)
        J2 = x**2 / (2D0*const_pi)**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (dexp(x)+2.3D0)
        J3 = 3D0*x**4 / (2D0*const_pi)**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (0.2D0 + x**2)
     .       / (dexp(x)+1D0+0.1D0*x**2)
        b0 =  3.4296D0
!             Possible error. I(0,E) should be same as he3_coll_int_ht
        b1 = -3.2148D0
        b2 =  2.375D0
        K0 = 9D0/8D0/(2D0*const_pi)**0.5D0
     .       / (1D0 + x)**0.5D0 !!! Was: (0.5D0 + x)**0.5D0
     .       / (dexp(x) + b0 + b1*x + b2*x**2)
!           Typo in the paper. This can be checked using Einzel-1978 f.80 high-temp
!           limit, I = 1 + (gap/ttc)^(A + B (Ep/T)^2) and calculating A and B. - slazav
        K1 = - 5D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**2 / (1D0+x)**0.5D0
     .       * dexp(-x) / (const_pi**2 + x**2)
        K2 = 3D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**2 / (1D0+x)**0.5D0
     .       * dexp(-x) / (127D0/150D0 * const_pi**2 + x**2)
        K3 = - 15D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**4 / (1D0+x)**0.5D0  !!! Was: x**4 / (1D0+x**2)**0.5D0
     .       * dexp(-x) / (const_pi**2 + x**2)**2
!           Possible typo in the paper: Ki sould be ~exp(x)/x at low temp. - slazav
        I0 = J0 + K0 * (xi/ttc)**2
        I1 = J1 + K1 * (xi/ttc)**2
        I2 = J2 + K2 * (xi/ttc)**2
        I3 = J3 + K3 * (xi/ttc)**2
        he3_coll_int = ( I0 - g0*(I1+I2) + d0*I3 )
      end

! Collision integral for low temp (good for < 0.7Tc)
! Einzel, JLTP84 (1991), p.345
      function he3_coll_int_lt(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 xi, ttc, gap, x, g0,d0,w0
        x=xi/dsqrt(2D0*ttc*gap)
        w0 = (1D0 - 2D0/3D0*g0 + d0)
        he3_coll_int_lt =
     .    3D0/2D0/const_pi * gap/ttc * he3_yosida(ttc,gap, 0D0)
     .    * (w0 + ttc/gap*(0.75D0*(1D0 + x**2) * w0
     .                      - (1D0+2D0*x**2)*(g0/3D0+d0) ))
      end

! Collision integral for high temp (good above 0.95 Tc)
! Einzel, JLTP84 (1991), p.345
! Einzel, JLTP32 (1978), f.80 - first (gap/ttc)^2 term
      function he3_coll_int_ht(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 xi, ttc, gap, x, g0,d0,w0
        he3_coll_int_ht = 1D0
     .     + (xi**2 + gap**2)/(ttc*const_pi)**2
     .     - (gap/ttc/const_pi)**2 *
     .          (6D0*dlog(2D0) + const_pi**2/18D0 
     .             + 3.3D0*(xi**2 + gap**2)/(ttc*const_pi)**2)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Quasiparticle lifetime at fermi level
! Einzel JLTP84 (1991) p.344
      function he3_tau0(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, g0, d0, tn
        if (ttc.lt.0D0) then
          he3_tau0 = NaN
          return
        endif
        if (ttc.gt.1D0) then
          he3_tau0=he3_tau_n0(ttc,p)
          return
        endif
        g0  = he3_scatt_g0(p)
        d0  = he3_scatt_d0(p)
        gap = he3_trivgap(ttc, p)
        tn  = he3_tau_n0(ttc, p)
        he3_tau0 = tn / he3_coll_int(0D0, ttc, gap, g0, d0);
      end

! Quasiparticle lifetime at low temp limit (no Ek dep)
! Einzel-1978 f.79
      function he3_tau0lt(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, g0, d0, tn
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_tau0 = NaN
          return
        endif
        gap = he3_trivgap(ttc, p)
        he3_tau0lt = he3_tau_n0(ttc, p) * dsqrt(2D0*const_pi) / 3D0
     .    * (ttc/gap)**1.5D0 * dexp(gap/ttc) / he3_scatt_w0(p)
      end

! Integrand for tau_av calculations
! Integration is similar to Y0 calculation in he3_gap.f
      function he3_tau_av_int(x, args)
        implicit none
        include 'he3.fh'
        real*8 x, args(4), ttc, gap, g0, d0, xi, Ek, C
        real*8 he3_tau_av_int

        ttc=args(1)
        gap=args(2)
        g0=args(3)
        d0=args(4)

        C=3D0 ! see tests/plot_tauav_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        he3_tau_av_int = he3_coll_int(xi,ttc, gap, g0, d0)
     .   / (dcosh(Ek/(2D0*ttc)))**2 / 2D0/ttc
     .   * C/(1D0-x**2)
      end

! Averaged quasiparticle lifetime
! Einzel JLTP84 (1991) p.345
      function he3_tau_av(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, g0, d0, Y0, tn, args(4)
        real*8 dx, xp, xm
        real*8 he3_tau_av_int
        external he3_tau_av_int
        if (ttc.lt.0D0) then
          he3_tau_av = NaN
          return
        endif
        if (ttc.gt.1D0) then
          he3_tau_av=he3_tau_n_av(ttc,p)
          return
        endif
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        tn  = he3_tau_n0(ttc, p)
        args = (/ttc, gap, he3_scatt_g0(p), he3_scatt_d0(p)/)
        he3_tau_av = Y0 * tn
     .   / math_dint(he3_tau_av_int, 0D0, 1D0, 100, args)
      end

! Integrands for he3_fpath calculation
! Integration is similar to Y0 calculation in he3_gap.f
      function he3_fpath_int1(x, args)
        implicit none
        include 'he3.fh'
        real*8 x,args(4),ttc, gap, g0, d0, xi, ek, I, C
        real*8 he3_fpath_int1
        ttc=args(1)
        gap=args(2)
        g0=args(3)
        d0=args(4)
        C=1D0 ! see tests/plot_tauav_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        I = he3_coll_int(xi,ttc, gap, g0, d0)
        he3_fpath_int1 = 1D0/I**2  ! (t/tN)^2
     .   * (xi/Ek)**2
     .   / (1D0 + dexp(Ek/ttc)) ! Fermi function
     .   * C/(1D0-x**2)
      end
      function he3_fpath_int2(x, args)
        implicit none
        include 'he3.fh'
        real*8 x,args(4),ttc, gap, xi, ek, C
        real*8 he3_fpath_int2
        ttc=args(1)
        gap=args(2)
        C=1D0 ! see tests/plot_tauav_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        he3_fpath_int2 =
     .   1D0 / (1D0 + dexp(Ek/ttc)) ! Fermi function
     .   * C/(1D0-x**2)
      end

! Mean free path of Bogoliubov quasiparticles
! Einzel JLTP32 (1978) f.84
      function he3_fpath(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, tn, vf, Ife, args(4)
        real*8 dx, xp, xm
        real*8 he3_fpath_int1, he3_fpath_int2
        external he3_fpath_int1, he3_fpath_int2
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_fpath=NaN
          return
        endif
        gap = he3_trivgap(ttc, p)
        tn  = he3_tau_n0(ttc, p)
        vf  = he3_vf(p)
        Ife = ttc * dlog(1D0+dexp(-gap/ttc)) ! integral of fermi function
        args = (/ttc,gap,he3_scatt_g0(p),he3_scatt_d0(p)/)
        he3_fpath = vf * tn
     .    *dsqrt(math_dint(he3_fpath_int1, 0D0, 1D0, 1000, args)/
     .           math_dint(he3_fpath_int2, 0D0, 1D0, 1000, args))
      end

! Spin diffusion perpendicular transport time, s
! Einzel JLTP84 (1991) f.90,96
      function he3_tau_dperp(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, l1a, gap, Y0, Yp
        l1a = he3_scatt_l1a(p)
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        Yp  = he3_yosida_perp(ttc, gap)
        he3_tau_dperp = he3_tau_av(ttc,p)
     .   /(1D0 - l1a*Yp/Y0)
      end

! Spin diffusion parallel transport time, s
! Einzel JLTP84 (1991) f.90,96
      function he3_tau_dpar(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, l1a, gap, Y0, Yp
        l1a = he3_scatt_l1a(p)
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        Yp  = he3_yosida(ttc, gap, 2D0)
        he3_tau_dpar = he3_tau_av(ttc,p)
     .   /(1D0 - l1a*Yp/Y0)
      end

! Hydrodynamic spin diffusion D_perp, cm2/s
! Einzel JLTP84 (1991) f.102
      function he3_diff_hperp_zz(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, tau, f0a, Y0, chi0, vf
        gap  = he3_trivgap(ttc, p)
        tau  = he3_tau_dperp(ttc,p)
        vf   = he3_vf(p)
        f0a  = he3_f0a(p)
        Y0   = he3_yosida(ttc, gap, 0D0)
        chi0 = (2D0+Y0)/(3D0+f0a*(2D0+Y0))
        he3_diff_hperp_zz = vf**2*tau/3D0/chi0
     .                    * he3_yosida_perp(ttc,gap)
      end

! Hydrodynamic spin diffusion D_par, cm2/s
! Einzel JLTP84 (1991) f.102
      function he3_diff_hpar_zz(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, tau, f0a, Y0, chi0, vf
        gap  = he3_trivgap(ttc, p)
        tau  = he3_tau_dpar(ttc,p)
        vf   = he3_vf(p)
        f0a  = he3_f0a(p)
        Y0   = he3_yosida(ttc, gap, 0D0)
        chi0 = (2D0+Y0)/(3D0+f0a*(2D0+Y0))
        he3_diff_hpar_zz = vf**2*tau/3D0/chi0
     .                 * he3_yosida_par(ttc,gap)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Integrand for spin diffusion calculation
! Integration is similar to Y0 calculation in he3_gap.f
! 2D integration needed: x [0:1] and th angle [0:pi]
! Bunkov, Dmitriev, Markelov, Mukharsky, PRL65 (1990) f.3
! Einzel, JLTP84 (1991) f.108 - typo in the formula: Szz- should be Sxx-
! Markelov, Mukharsky, Ph.B178 (1992) - most general result
! arguments:
!   ttc
!   gap
!   o0     -- Larmor freq (rad/s)
!   lambda -- exchange coupling strength
!   td     -- Spin diffusion transport time, s
!   type   1: D_perp_zz
!          2: D_perp_xx

      function he3_diff_int(x, th, args)
        implicit none
        real*8 x,th, args(6), ttc, gap, o0, lambda, td
        real*8 C, xi, Ek, kz, kp, phi, u, v, o1, Sp2, Sm2
        integer type
        complex*16 t,s
        complex*16 he3_diff_int

        ttc=args(1)
        gap=args(2)
        o0     = args(3)
        lambda = args(4)
        td     = args(5)
        type   = nint(args(6))

        C=2D0
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        phi = (dcosh(Ek/(2D0*ttc)))**(-2) / 2D0 / ttc
        u = xi/Ek
        v = gap/Ek

        kz=dsin(th)
        kp=dcos(th)
        o1 = o0*(1+lambda)

        he3_diff_int = (0D0,0D0)

        Sm2 = 1D0 - kp**2/2 *(1-u**2)
        Sp2 = u**2 + (1-u**2)*kz**2
        t = dcmplx(td, 0D0) / dcmplx(1D0, -o0*td)
        s = dcmplx(o1, 0D0) * t

        if (type.eq.1) then ! D_perp_xx
          he3_diff_int = t * dcmplx(0.5D0*kp * kp**2, 0D0)
     .      * (Sm2 - Sp2*s * (0D0,1D0)) / (1D0 + Sp2*s**2)
     .      * dcmplx(phi * C/(1D0-x**2), 0D0)
          return
        endif

        if (type.eq.2) then ! D_perp_zz
          he3_diff_int = t * dcmplx(kp *kz**2, 0D0)
     .      * (Sm2 - Sp2*s * (0D0,1D0)) / (1D0 + Sp2*s**2)
     .      * dcmplx(phi * C/(1D0-x**2), 0D0)
          return
        endif

        if (type.eq.3) then ! D_par_xx
          he3_diff_int = dcmplx( td * 0.5D0*kp * kp**2
     .      * (Sm2 + u**2 * (o1*td)**2) / (1 + Sp2*(o1*td)**2)
     .      * phi * C/(1D0-x**2), 0D0)
          return
        endif

        if (type.eq.4) then ! D_par_zz
          he3_diff_int = dcmplx( td * kp * kz**2
     .      * (Sm2 + u**2 * (o1*td)**2) / (1D0 + Sp2*(o1*td)**2)
     .      * phi * C/(1D0-x**2), 0D0)
          return
        endif

      end

! Spin diffusion coefficient D_perp, cm2/s
      function he3_diff_all(ttc, p, nu0, type)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0, type
        real*8 args(6), gap, Vf, chi0, Y0, f0a
        real*8 o0, lambda, td
        complex*16 he3_diff_int
        external he3_diff_int
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_diff_all=NaN
          return
        endif
        gap = he3_trivgap(ttc, p)
        Vf  = he3_vf(p)
        f0a = he3_f0a(p)
        Y0  = he3_yosida(ttc, gap, 0D0);
        chi0    = (2D0 + Y0) / (3D0 + f0a*(2D0 + Y0))
        o0      = nu0*2D0*const_pi
        lambda  = -f0a*chi0  ! Einzel-1991 p.349
        td      = he3_tau_dperp(ttc, p)

        if (type.lt.10D0) then
          args = (/ttc, gap, o0, lambda, td, type/)
          he3_diff_all = dreal(math_cint2d(he3_diff_int,
     .          0D0, 1D0, 100, 0D0, const_pi/2D0, 100, args))
     .      * Vf**2 / chi0
        else
          args = (/ttc, gap, o0, lambda, td, type-10D0/)
          he3_diff_all = dimag(math_cint2d(he3_diff_int,
     .          0D0, 1D0, 100, 0D0, const_pi/2D0, 100, args))
     .      * Vf**2 / chi0
        endif
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spin diffusion coefficient D_perp_xx, cm2/s
      function he3_diff_perp_xx(ttc, p, nu0)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        if (ttc.ge.1D0) then
          he3_diff_perp_xx = he3_diffn_perp(ttc,p,nu0)
        else
          he3_diff_perp_xx = he3_diff_all(ttc, p, nu0, 1D0)
        endif
      end

! Spin diffusion coefficient D_perp_xx_im, cm2/s
      function he3_diff_perp_xx_im(ttc, p, nu0)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        if (ttc.ge.1D0) then
          he3_diff_perp_xx_im = NaN
        else
          he3_diff_perp_xx_im = he3_diff_all(ttc, p, nu0, 11D0)
        endif
      end

! Spin diffusion coefficient D_perp_zz, cm2/s
      function he3_diff_perp_zz(ttc, p, nu0)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        if (ttc.ge.1D0) then
          he3_diff_perp_zz = he3_diffn_perp(ttc,p, nu0)
        else
          he3_diff_perp_zz = he3_diff_all(ttc, p, nu0, 2D0)
        endif
      end

! Spin diffusion coefficient D_perp_zz_im, cm2/s
      function he3_diff_perp_zz_im(ttc, p, nu0)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        if (ttc.ge.1D0) then
          he3_diff_perp_zz_im = NaN
        else
          he3_diff_perp_zz_im = he3_diff_all(ttc, p, nu0, 12D0)
        endif
      end

! Spin diffusion coefficient D_par_xx, cm2/s
      function he3_diff_par_xx(ttc, p, nu0)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        if (ttc.ge.1D0) then
          he3_diff_par_xx = he3_diffn_hydr(ttc,p)
        else
          he3_diff_par_xx = he3_diff_all(ttc, p, nu0, 3D0)
        endif
      end

! Spin diffusion coefficient D_perp_zz, cm2/s
      function he3_diff_par_zz(ttc, p, nu0)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        if (ttc.ge.1D0) then
          he3_diff_par_zz = he3_diffn_hydr(ttc,p)
        else
          he3_diff_par_zz = he3_diff_all(ttc, p, nu0, 4D0)
        endif
      end

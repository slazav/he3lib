!HH> He3-B transport properties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Collision integral for Bogoliubov quasiparticles

!> Collision integral in Einzel approximation
!> Einzel, Wolfle, Hirschfeld, JLTP80 (1990), Appendix, p.66
!> + my small fixes
      function he3_coll_int(xi,ttc, gap, g0, d0) !F>
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
        J1 = x**2 / const_2pi**0.5D0
     .       * (0.5D0 + x)**1.5D0
     .       / (1.3D0 + x**2) * dexp(-x)
        J2 = x**2 / const_2pi**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (dexp(x)+2.3D0)
        J3 = 3D0*x**4 / const_2pi**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (0.2D0 + x**2)
     .       / (dexp(x)+1D0+0.1D0*x**2)
        b0 =  3.4296D0
!             Possible error. I(0,E) should be same as he3_coll_int_ht
        b1 = -3.2148D0
        b2 =  2.375D0
        K0 = 9D0/8D0/const_2pi**0.5D0
     .       / (1D0 + x)**0.5D0 !!! Was: (0.5D0 + x)**0.5D0
     .       / (dexp(x) + b0 + b1*x + b2*x**2)
!           Typo in the paper. This can be checked using Einzel-1978 f.80 high-temp
!           limit, I = 1 + (gap/ttc)^(A + B (Ep/T)^2) and calculating A and B. - slazav
        K1 = - 5D0/8D0/const_2pi**0.5D0
     .       * x**2 / (1D0+x)**0.5D0
     .       * dexp(-x) / (const_pi**2 + x**2)
        K2 = 3D0/8D0/const_2pi**0.5D0
     .       * x**2 / (1D0+x)**0.5D0
     .       * dexp(-x) / (127D0/150D0 * const_pi**2 + x**2)
        K3 = - 15D0/8D0/const_2pi**0.5D0
     .       * x**4 / (1D0+x)**0.5D0  !!! Was: x**4 / (1D0+x**2)**0.5D0
     .       * dexp(-x) / (const_pi**2 + x**2)**2
!           Possible typo in the paper: Ki should be ~exp(x)/x at low temp. - slazav
        I0 = J0 + K0 * (xi/ttc)**2
        I1 = J1 + K1 * (xi/ttc)**2
        I2 = J2 + K2 * (xi/ttc)**2
        I3 = J3 + K3 * (xi/ttc)**2
        he3_coll_int = ( I0 - g0*(I1+I2) + d0*I3 )
      end

!> Collision integral for low temp (good for < 0.7Tc)
!> Einzel, JLTP84 (1991), p.345
      function he3_coll_int_lt(xi,ttc, gap, g0, d0) !F>
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

!> Collision integral for high temp (good above 0.95 Tc)
!> Einzel, JLTP84 (1991), p.345
!> Einzel, JLTP32 (1978), f.80 - first (gap/ttc)^2 term
      function he3_coll_int_ht(xi,ttc, gap, g0, d0) !F>
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

!H> Quasiparticle lifetime, mean free path

!> Quasiparticle lifetime at Fermi level
!> Einzel JLTP84 (1991) p.344
      function he3_tau0(ttc, p) !F>
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
        gap = he3_gap(ttc, p)
        tn  = he3_tau_n0(ttc, p)
        he3_tau0 = tn / he3_coll_int(0D0, ttc, gap, g0, d0);
      end

!> Quasiparticle lifetime at low temp limit (no Ek dep)
!> Einzel-1978 f.79
      function he3_tau0lt(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, g0, d0, tn
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_tau0 = NaN
          return
        endif
        gap = he3_gap(ttc, p)
        he3_tau0lt = he3_tau_n0(ttc, p) * dsqrt(const_2pi) / 3D0
     .    * (ttc/gap)**1.5D0 * dexp(gap/ttc) / he3_scatt_w0(p)
      end

! Integrand for tau_av calculations
! Integration is similar to Y0 calculation in he3_gap.f
      function he3_tau_av_int(x)
        implicit none
        include 'he3.fh'
        real*8 x, he3_tau_av_int
        real*8 ttc, gap, g0, d0
        common /he3_tau_av_int_cb/ ttc, gap, g0, d0
        real*8 xi, Ek, C

        C=3D0 ! see tests/plot_tauav_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        he3_tau_av_int = he3_coll_int(xi,ttc, gap, g0, d0)
     .   / (dcosh(Ek/(2D0*ttc)))**2 / 2D0/ttc
     .   * C/(1D0-x**2)
      end

!> Averaged quasiparticle lifetime
!> Einzel JLTP84 (1991) p.345
      function he3_tau_av(ttc, p) !F>
        implicit none
        include 'he3.fh'
        include 'he3_math.fh'
        real*8 ttc, p
        real*8 he3_tau_av_int
        external he3_tau_av_int
        real*8 ttc1, gap, g0, d0, Y0, tn
        common /he3_tau_av_int_cb/ ttc1, gap, g0, d0

        if (ttc.lt.0D0) then
          he3_tau_av = NaN
          return
        endif
        if (ttc.gt.1D0) then
          he3_tau_av=he3_tau_n_av(ttc,p)
          return
        endif
        ttc1=ttc
        gap = he3_gap(ttc, p)
        g0  = he3_scatt_g0(p)
        d0  = he3_scatt_d0(p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        tn  = he3_tau_n0(ttc, p)
        he3_tau_av = Y0 * tn
     .   / math_dint(he3_tau_av_int, 0D0, 1D0, 100)
      end

! Integrands for he3_fpath calculation
! Integration is similar to Y0 calculation in he3_gap.f
      function he3_fpath_int1(x)
        implicit none
        include 'he3.fh'
        real*8 x, he3_fpath_int1
        real*8 ttc, gap, g0, d0
        common /he3_fpath_int1_cb/ ttc, gap, g0, d0
        real*8 xi, ek, I, C

        C=1D0 ! see tests/plot_tauav_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        I = he3_coll_int(xi,ttc, gap, g0, d0)
        he3_fpath_int1 = 1D0/I**2  ! (t/tN)^2
     .   * (xi/Ek)**2
     .   / (1D0 + dexp(Ek/ttc)) ! Fermi function
     .   * C/(1D0-x**2)
      end
      function he3_fpath_int2(x)
        implicit none
        include 'he3.fh'
        real*8 x, he3_fpath_int2
        real*8 ttc, gap
        common /he3_fpath_int2_cb/ ttc, gap
        real*8 xi, ek, C

        C=1D0 ! see tests/plot_tauav_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        he3_fpath_int2 =
     .   1D0 / (1D0 + dexp(Ek/ttc)) ! Fermi function
     .   * C/(1D0-x**2)
      end

!> Mean free path of Bogoliubov quasiparticles [cm]
!> Einzel JLTP32 (1978) f.84
      function he3_fpath(ttc, p) !F>
        implicit none
        include 'he3.fh'
        include 'he3_math.fh'
        real*8 ttc, p, tn, vf
        real*8 he3_fpath_int1, he3_fpath_int2
        external he3_fpath_int1, he3_fpath_int2
        real*8 ttc1, gap1, g0, d0, ttc2, gap2
        common /he3_fpath_int1_cb/ ttc1, gap1, g0, d0
        common /he3_fpath_int2_cb/ ttc2, gap2

        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_fpath=NaN
          return
        endif
        ttc1=ttc
        ttc2=ttc
        gap1 = he3_gap(ttc, p)
        gap2 = gap1
        g0  = he3_scatt_g0(p)
        d0  = he3_scatt_d0(p)
        tn  = he3_tau_n0(ttc, p)
        vf  = he3_vf(p)
        he3_fpath = vf * tn
     .    *dsqrt(math_dint(he3_fpath_int1, 0D0, 1D0, 1000)/
     .           math_dint(he3_fpath_int2, 0D0, 1D0, 1000))
      end

!> RMS group velocity of Bogoliubov quasiparticles [cm/s]
!> Einzel JLTP32 (1990) f.28 and below
      function he3_rmsv(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, Y0,Y2

        gap = he3_gap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        Y2  = he3_yosida(ttc, gap, 2D0)
        he3_rmsv = he3_vf(p) * dsqrt(Y2/Y0)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Viscosity

!> Viscous free path of Bogoliubov quasiparticles [cm]
!> Einzel 1990 Eq.26
      function he3_visc_fpath(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, Y0,Y2, l2

        ! qubic fit of l2 from Einzel-1990, table.1:
        l2 = 5D-6*p**3 - 4D-4*p**2 + 9.5D-3*p + 0.68D0

        gap = he3_gap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        Y2  = he3_yosida(ttc, gap, 2D0)

        he3_visc_fpath = he3_fpath(ttc, p) / (1D0 - l2*Y2/Y0)
      end

!> Hydrodinamic (freq=0) viscosity of Bogoliubov quasiparticles [g/cm/s]
!> Einzel 1990 Eq.28
      function he3_hvisc(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, Y0,Y2, l2

        gap = he3_gap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        Y2  = he3_yosida(ttc, gap, 2D0)

        ! Note: in Eq.28 rho0n = rho*Y is used
        he3_hvisc = 0.2D0 * he3_rho(p)
     .   * he3_vf(p) * dsqrt(Y2*Y0)
     .   * he3_visc_fpath(ttc,p)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Spin diffusion

!> Spin diffusion perpendicular transport time, s
!> Einzel JLTP84 (1991) f.90,96
      function he3_tau_dperp(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, l1a, gap, Y0, Yp
        l1a = he3_scatt_l1a(p)
        gap = he3_gap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        Yp  = he3_yosida_perp(ttc, gap)
        he3_tau_dperp = he3_tau_av(ttc,p)
     .   /(1D0 - l1a*Yp/Y0)
      end

!> Spin diffusion parallel transport time, s
!> Einzel JLTP84 (1991) f.90,96
      function he3_tau_dpar(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, l1a, gap, Y0, Yp
        l1a = he3_scatt_l1a(p)
        gap = he3_gap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        Yp  = he3_yosida(ttc, gap, 2D0)
        he3_tau_dpar = he3_tau_av(ttc,p)
     .   /(1D0 - l1a*Yp/Y0)
      end

!> Hydrodynamic spin diffusion D_perp, cm2/s
!> Einzel JLTP84 (1991) f.102
      function he3_diff_hperp_zz(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, tau, f0a, Y0, chi0, vf
        gap  = he3_gap(ttc, p)
        tau  = he3_tau_dperp(ttc,p)
        vf   = he3_vf(p)
        f0a  = he3_f0a(p)
        Y0   = he3_yosida(ttc, gap, 0D0)
        chi0 = (2D0+Y0)/(3D0+f0a*(2D0+Y0))
        he3_diff_hperp_zz = vf**2*tau/3D0/chi0
     .                    * he3_yosida_perp(ttc,gap)
      end

!> Hydrodynamic spin diffusion D_par, cm2/s
!> Einzel JLTP84 (1991) f.102
      function he3_diff_hpar_zz(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, tau, f0a, Y0, chi0, vf
        gap  = he3_gap(ttc, p)
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

      function he3_diff_int(x, kz)
        implicit none
        real*8 x,kz, he3_diff_int

        real*8 ttc, gap, o0, lambda, td
        integer itype
        common /he3_diff_int_cb/ ttc, gap, o0, lambda, td, itype

        real*8 C, xi, Ek, phi, uu, vv, o1, Sp2, Sm2
        complex*16 i,e,t,s,res, cSp2,cSm2

        C=3.5D0*ttc ! see plot_sdiff_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        phi = (dcosh(Ek/(2D0*ttc)))**(-2) / 2D0 / ttc
        uu = (xi/Ek)**2
        vv = (gap/Ek)**2

        o1 = o0*(1D0+lambda)

        Sm2 = 1D0 - (1D0-kz**2)/2D0 * vv
        Sp2 = uu + vv*kz**2

        i = dcmplx(0D0,1D0)
        e = dcmplx(1D0,0D0)
        t = dcmplx(td, 0D0) / dcmplx(1D0, -o0*td)
        s = dcmplx(o1, 0D0) * t

        cSm2 = dcmplx(Sm2,0D0)
        cSp2 = dcmplx(Sp2,0D0)

        res=(0D0,0D0)
        if (itype.eq.1.or.itype.eq.11) then ! D_perp_xx
          res = t * dcmplx((1D0-kz**2)/2D0, 0D0)
     .      * (cSm2 - cSp2*s*i) / (e + cSp2*s**2)
        elseif (itype.eq.2.or.itype.eq.12) then ! D_perp_zz
          res = t * dcmplx(kz**2, 0D0)
     .      * (cSm2 - cSp2*s*i) / (e + cSp2*s**2)
        elseif (itype.eq.3.or.itype.eq.13) then ! D_par_xx
          res = dcmplx( td * (1D0-kz**2)/2D0
     .      * (Sm2 + uu*(o1*td)**2)/(1D0 + Sp2*(o1*td)**2), 0D0)
        elseif (itype.eq.4.or.itype.eq.14) then ! D_par_zz
          res = dcmplx( td * kz**2
     .      * (Sm2 + uu*(o1*td)**2)/(1D0 + Sp2*(o1*td)**2), 0D0)
        endif

        if (itype.lt.10) then
          he3_diff_int = dreal(res)
        else
          he3_diff_int = dimag(res)
        endif
        he3_diff_int = he3_diff_int * phi * C/(1D0-x**2)
      end

! same, but integrated over kz
! see plot_sdiff_int1.m
      function he3_diff_int_i(x)
        implicit none
        real*8 x, he3_diff_int_i

        real*8 ttc, gap, o0, lambda, td
        integer itype
        common /he3_diff_int_cb/ ttc, gap, o0, lambda, td, itype

        real*8 C, xi, Ek, phi, o1
        complex*16 i,e,h, t,s,res, uu, vv, ot2
        complex*16 AA, BB, CC, DD, AC,CD,D1,T1, I1,I2

        C=4.0D0 ! see plot_sdiff_int.m
        xi = datanh(x)*C
        Ek=dsqrt(xi**2 + gap**2)
        phi = (dcosh(Ek/(2D0*ttc)))**(-2) / 2D0 / ttc
        uu = dcmplx((xi/Ek)**2, 0D0)
        vv = dcmplx((gap/Ek)**2, 0D0)

        o1 = o0*(1D0+lambda)
        ot2=dcmplx((o1*td)**2, 0D0)

        i = dcmplx(0D0,1D0)
        e = dcmplx(1D0,0D0)
        h = dcmplx(0.5D0,0D0)
        t = dcmplx(td, 0D0) / dcmplx(1D0, -o0*td)
        s = dcmplx(o1, 0D0) * t

        if (itype.eq.1.or.itype.eq.11.or.
     .      itype.eq.2.or.itype.eq.12) then ! D_perp
          AA = (h - i*s)*vv
          BB = h*(e+uu) - i*s*uu
          CC = s**2 * vv
          DD = e + s**2 * uu
        else ! D_parallel
          AA = h*vv
          BB = h*(e+uu) + ot2 * uu
          CC = ot2 * vv
          DD = e + ot2 * uu
        endif

        if (abs(CC)<1D-4) then ! close to hydrodynamic OR x close to 1
          I1 = (AA/(3D0,0D0) + BB
     .       - AA*CC/DD/(5D0,0D0) - BB*CC/DD/(3D0,0D0))/DD
          I2 = (AA/(5D0,0D0) + BB/(3D0,0D0)
     .       - AA*CC/DD/(7D0,0D0) - BB*CC/DD/(5D0,0D0))/DD
        else
          AC = AA/CC
          CD = CC/DD
          D1 = BB/AA - DD/CC
          ! atan x = 1/2 i (log(1-ix)-log(1+ix))
          T1 = cdsqrt(CD) *
     .       h*i*(cdlog(e-i*cdsqrt(CD)) - cdlog(e+i*cdsqrt(CD)))
          I1 = AC + AC*D1*T1
          I2 = AC/(3D0,0D0) + AC*D1*(e-T1/CD)
        endif

        res=(0D0,0D0)
        if (itype.eq.1.or.itype.eq.11) then ! D_perp_xx
          res = t * (I1-I2)/(2D0,0D0)
        elseif (itype.eq.2.or.itype.eq.12) then ! D_perp_zz
          res = t * I2
        elseif (itype.eq.3.or.itype.eq.13) then ! D_par_xx
          res = dcmplx(td,0D0) * (I1-I2)/(2D0,0D0)
        elseif (itype.eq.4.or.itype.eq.14) then ! D_par_zz
          res = dcmplx(td,0D0) * I2
        endif

        if (itype.lt.10) then
          he3_diff_int_i = dreal(res)
        else
          he3_diff_int_i = dimag(res)
        endif
        he3_diff_int_i = he3_diff_int_i * phi * C/(1D0-x**2)
      end

! diffusion coefficient D, any component depending on type [cm2/s]
      function he3_diff_all(ttc, p, nu0, type)
        implicit none
        include 'he3.fh'
        include 'he3_math.fh'
        real*8 ttc, p, nu0, type
        real*8 Vf, chi0, Y0, f0a
        real*8 he3_diff_int_i, he3_diff_int
        external he3_diff_int_i, he3_diff_int, math_dint_gka
        real*8 he3_diff_all
        real*8 rerr,aerr

        real*8 ttc1,gap,o0, lambda, td
        integer itype
        common /he3_diff_int_cb/ ttc1, gap, o0, lambda, td, itype

        if (ttc.lt.0D0.or.ttc.gt.1D0.or.isnan(ttc)) then
          he3_diff_all=NaN
          return
        endif
        if (ttc.lt.1D-2) then ! less then 1D-70
          he3_diff_all=0D0
          return
        endif

        ttc1 = ttc
        gap = he3_gap(ttc, p)
        Vf  = he3_vf(p)
        f0a = he3_f0a(p)
        Y0  = he3_yosida(ttc, gap, 0D0);
        chi0    = (2D0 + Y0) / (3D0 + f0a*(2D0 + Y0))
        o0      = nu0*const_2pi
        lambda  = -f0a*chi0  ! Einzel-1991 p.349

        itype   = nint(type)

        if (itype.eq.1.or.itype.eq.11.or.
     .      itype.eq.2.or.itype.eq.12) then
          td = he3_tau_dperp(ttc, p)
        else
          td = he3_tau_dpar(ttc, p)
        endif

!       slow 2D integration
!        he3_diff_all = math_dint2d(he3_diff_int,
!     .    0D0, 1D0, 200, 0D0,1D0, 200)
!     .    * Vf**2 / chi0

!       fast 1D integration (kz is integrated analytically)
!        he3_diff_all = math_dint(he3_diff_int_i, 0D0, 1D0, 500)
!     .    * Vf**2 / chi0
        rerr=1D-8
        aerr=1D-40
        he3_diff_all = math_dint_gka(
     .                 he3_diff_int_i, 0D0, 1D0, 0D0, rerr)
     .                 * Vf**2 / chi0

!        he3_diff_all=0D0
!        call math_dint2d_ad(math_dint2d_ad, he3_diff_int,
!     .    0D0, 1D0, 0D0, const_pi/2D0, 0D0, 1D-5, he3_diff_all)
!        he3_diff_all = he3_diff_all * Vf**2 / chi0
!       write(*,*) 'int>', he3_diff_all
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Spin diffusion coefficient D_perp_xx [cm2/s]
      function he3_diff_perp_xx(ttc, p, nu0) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        real*8 he3_diff_all
        if (ttc.ge.1D0) then
          he3_diff_perp_xx = he3_diffn_perp(ttc,p,nu0)
        else
          he3_diff_perp_xx = he3_diff_all(ttc, p, nu0, 1D0)
        endif
      end

!> Spin diffusion coefficient D_perp_xx_im [cm2/s]
      function he3_diff_perp_xx_im(ttc, p, nu0) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        real*8 he3_diff_all
        if (ttc.ge.1D0) then
          he3_diff_perp_xx_im = NaN
        else
          he3_diff_perp_xx_im = he3_diff_all(ttc, p, nu0, 11D0)
        endif
      end

!> Spin diffusion coefficient D_perp_zz [cm2/s]
      function he3_diff_perp_zz(ttc, p, nu0) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        real*8 he3_diff_all
        if (ttc.ge.1D0) then
          he3_diff_perp_zz = he3_diffn_perp(ttc,p, nu0)
        else
          he3_diff_perp_zz = he3_diff_all(ttc, p, nu0, 2D0)
        endif
      end

!> Spin diffusion coefficient D_perp_zz_im [cm2/s]
      function he3_diff_perp_zz_im(ttc, p, nu0) !F>
        implicit none
        include 'he3.fh'
        real*8 he3_diff_all
        real*8 ttc, p, nu0
        if (ttc.ge.1D0) then
          he3_diff_perp_zz_im = NaN
        else
          he3_diff_perp_zz_im = he3_diff_all(ttc, p, nu0, 12D0)
        endif
      end

!> Spin diffusion coefficient D_par_xx [cm2/s]
      function he3_diff_par_xx(ttc, p, nu0) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        real*8 he3_diff_all
        if (ttc.ge.1D0) then
          he3_diff_par_xx = he3_diffn_hydr(ttc,p)
        else
          he3_diff_par_xx = he3_diff_all(ttc, p, nu0, 3D0)
        endif
      end

!> Spin diffusion coefficient D_perp_zz [cm2/s]
      function he3_diff_par_zz(ttc, p, nu0) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        real*8 he3_diff_all
        if (ttc.ge.1D0) then
          he3_diff_par_zz = he3_diffn_hydr(ttc,p)
        else
          he3_diff_par_zz = he3_diff_all(ttc, p, nu0, 4D0)
        endif
      end


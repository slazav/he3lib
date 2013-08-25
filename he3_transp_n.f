! Normal phase spin diffusion

! Crossections

! Ti, Si parameters (for using in <W> calculations)
! Einzel & Wolfle JLTP32 (1987) page 34
      subroutine he3_s0s1t0t1(P, S0,S1,T0,T1)
        implicit none
        include 'he3.fh'
        real*8 P
        real*8 f0s, f0a, f1s, f1a
        real*8 A0s, A0a, A1s, A1a
        real*8 S0,S1,T0,T1, W, Wa
        f0s = he3_f0s(P)
        f0a = he3_f0a(P)
        f1s = he3_f1s(P)
        f1a = he3_f1a(P)

        A0s = f0s/(1D0+f0s)
        A0a = f0a/(1D0+f0a)
        A1s = f1s/(1D0+f1s/3D0)
        A1a = f1a/(1D0+f1a/3D0)

        S0 = A0s - 3D0*A0a
        S1 = A1s - 3D0*A1a
        T0 = A0s + A0a
        T1 = A1s + A1a
      end

! Scattering crossection
! Einzel & Wolfle JLTP32 (1978) f.82
      function he3_crsect_w(P)
        implicit none
        include 'he3.fh'
        real*8 P, S0,S1,T0,T1
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        he3_crsect_w = const_pi/2D0 *
     .    (S0**2 - 2D0/3D0*S0*S1 + 7D0/15D0*S1**2
     .     + 1.5D0*T0**2 - T0*T1 + 7D0/10D0*T1**2)
      end

! Einzel & Wolfle JLTP32 (1978) f.82
      function he3_crsect_wi(P)
        implicit none
        include 'he3.fh'
        real*8 P, S0,S1,T0,T1
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        he3_crsect_wi = const_pi/2D0 *
     .    ((S0**2)/3D0 + 2D0/15D0*S0*S1 - 29D0/105D0*S1**2
     .     + 2D0/3D0*S0*T0 - 2D0/5D0*S0*T1
     .     + 5D0/3D0*(25D0-36D0*dlog(2D0))*T0**2
     .     + (84D0-120D0*dlog(2D0))*T0*T1
     .     + 5D0/21D0*(173D0-252D0*dlog(2D0))*T1**2)
      end
      function he3_crsect_wd(P)
        implicit none
        include 'he3.fh'
        real*8 P, S0,S1,T0,T1
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        he3_crsect_wd = const_pi/2D0 *
     .    (7D0/15D0*S0**2 - 18D0/35D0*S0*S1 + 107D0/315D0*S1**2
     .     + 8D0/15D0*S0*T0 - 8D0/105D0*(S0*T1 + S1*T0)
     .     + 8D0/63D0*S1*T1 + 29D0/30D0*T0**2
     .     - 19D0/35D0*T0*T1 + 33D0/70D0*T1**2)
      end
! Einzel & Wolfle JLTP32 (1978) f.74
! code from Samuli
      function he3_crsect_wl(P)
        implicit none
        include 'he3.fh'
        real*8 P, S0,S1,T0,T1
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        he3_crsect_wl = const_pi * 1D0/420D0 *
     .    (- 70D0*S0**2 - 54D0*S1**2 + 175D0*T0**2
     .     + 28D0*S0*(7D0*S1 + 10D0*T0 - 6D0*T1)
     .     - 42D0*T0*T1 + 71D0*T1**2
     .     + 8D0*S1*(-21D0*T0 + 19D0*T1))
      end

! Scattering factors

! Einzel & Wolfle JLTP32 (1978) f.74
! As noticed in Einzel-91 p.350, \lambda_1^a can
! be set to 0.
      function He3_scatt_l1a(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_scatt_l1a =
     .    he3_crsect_wl(P) / he3_crsect_w(P)
      end
! Einzel & Wolfle JLTP32 (1978) f.66
      function He3_scatt_g0(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_scatt_g0 =
     .    he3_crsect_wi(P) / he3_crsect_w(P)
      end
! Einzel & Wolfle JLTP32 (1978) f.67
      function He3_scatt_d0(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_scatt_d0 =
     .    he3_crsect_wd(P) / he3_crsect_w(P)
      end
      function He3_scatt_w0(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_scatt_w0 = 1D0
     .    - 2D0/3D0 * He3_scatt_g0(P)
     .    + He3_scatt_d0(P)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Normal state quasiparticle lifetime at the Fermi level
! \tau_N(0,T)
! Einzel JLTP32 (1978) p.28,34
! Einzel JLTP84 (1991) f.4
      function He3_tau_n0(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        real*8 S0,S1,T0,T1, W
        W = he3_crsect_w(p)
        He3_tau_n0 = 32D0 *
     .    he3_tfeff(P)*const_hbar/const_pi**2
     .    / W / const_kb / (1D-3*ttc*he3_tc(P))**2
      end

! Thermal average of normal state quasiparticle lifetime
! \bar\tau_N(T)
! Einzel JLTP84 (1991) f.5
      function He3_tau_n_av(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        He3_tau_n_av = 0.75D0 * He3_tau_n0(ttc, p)
      end

! Spin diffusion transport time for a normal Fermi-liquid, s
! Einzel JLTP84 (1991) p.328
      function He3_tau_nd(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        He3_tau_nd = 0.75D0 * He3_tau_n0(ttc, p)
     .    / (1D0-he3_scatt_l1a(p))
      end

! Hydrodynamic spin diffusion in normal liquid, cm2/s
! Einzel JLTP84 (1991) f.23
      function he3_diffn_hydr(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, f0a, tau, vf
        tau  = he3_tau_nd(ttc,p)
        vf   = he3_vf(p)
        f0a  = he3_f0a(p)
        he3_diffn_hydr = vf**2 / 3D0 * (1D0+f0a) * tau
      end

! Spin diffusion Dperp in normal liquid, cm2/s
! Einzel JLTP84 (1991) f.22, Bunkov PRL65
      function he3_diffn_perp(ttc, p, nu0)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        real*8 vf, f0a, tau, oe
        f0a  = he3_f0a(p)
        vf   = he3_vf(p)
        tau  = he3_tau_nd(ttc,p)
        oe  = -f0a/(1D0+f0a) * nu0*2D0*const_pi
        he3_diffn_perp = vf**2 / 3D0 * (1D0+f0a)
     .    * tau /(1D0 + (tau * oe)**2)
      end

!HH> Normal He3 transport properties

!H> Crossections

! s-p approximation, only l=0,1 fermi-liquid parameters are non-zero.
! then (see Einzel-1978 f.81):
!  As(th, phi) = S0 + S1 cos(th)
!  At(th, phi) = (T0 + T1 cos(th))cos(ph)
!  A2s(th, phi) = S0 + S1 cos(th2)
!  A2t(th, phi) = (T0 + T1 cos(th2))cos(ph2)
!
! S and T are given by following subroutine:

! Ti, Si parameters (for using in <W> calculations)
! Einzel & Wolfle JLTP32 (1978) page 34
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

! Average over angles:
!  <...>  = int(0, pi) sin(th)/2cos(th/2) dth  int (0..2pi) dphi/2pi ...
!  <1> = 2
!  <cos(th)> = -2/3
!  <cos^2(th)> = 14/15
!  <cos^3(th)> = -18/35
!  <cos^2(ph)> = 1
!  <As> = 2S0 - S1 2/3, <At> = 0
!
! W(th, phi) = pi/4 (3 At^2 + As^2)
! = pi/4 <
!  + S0^2
!  + 2 S0 S1 cos(th)
!  + S1^2 cos^2(th)
!  + 3 T0^2 cos^2(ph)
!  + 6 T0 T1 cos(th)cos^2(ph)
!  + 3 T1^2 cos^2(th)cos^2(ph)
! >

! <W> =
! = pi/2 <
!  + 1    S0^2
!  - 2/3  S0 S1
!  + 7/15 S1^2
!  + 3/2  T0^2
!  - 1    T0 T1
!  + 7/10 T1^2
! >

! <W*cos(th)> =
! = pi/2 <
!  - 1/3   S0^2
!  + 14/15 S0 S1
!  - 9/35  S1^2
!  - 1/2   T0^2
!  + 7/5   T0 T1
!  - 27/70 T1^2
! >

! <W*sin^4(th/2)sin^2(phi)> =
! = pi/2 <
!  + 4/15    S0^2
!  - 8/21    S0 S1
!  + 52/315  S1^2
!  + 1/5     T0^2
!  - 2/7     T0 T1
!  + 13/105  T1^2
! >


!> Scattering crossection <W>, Einzel & Wolfle JLTP32 (1978) f.82
      function he3_crsect_w(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P, S0,S1,T0,T1
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        he3_crsect_w = const_pi/2D0 *
     .    (S0**2 - 2D0/3D0*S0*S1 + 7D0/15D0*S1**2
     .     + 1.5D0*T0**2 - T0*T1 + 7D0/10D0*T1**2)
      end

!> Scattering crossection <Wi>, Einzel & Wolfle JLTP32 (1978) f.82
      function he3_crsect_wi(P) !F>
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

!> Scattering crossection <Wd>, Einzel & Wolfle JLTP32 (1978) f.82
      function he3_crsect_wd(P) !F>
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

!> <p>Scattering parameters
!> $\lambda_n^+$ ($\lambda_n$, $\lambda_n^s$), $\lambda_n^-$ ($\lambda_n^a$),
!> $\delta_n^+$, $\delta_n^-$,
!> $\gamma_n^+$, $\gamma_n^-$,
!> <br>$\lambda_0^+ = \lambda_1^+ = 1$, $\lambda_0^- = 3$,
!> $\delta_0^+ = \delta_1^+ = \delta_0^-/3 = \delta_0$
!> <br>See Sykes-1970 f.26-29, Einzel-1978 f66,67,71,74, Einzel-1984 f.24


!> Scattering parameter $\lambda_1^-$ ($\lambda_1^a$)
!> l1a = 1 + 2 &lt;W*cos(th)&gt;/&lt;W&gt;
      function he3_scatt_l1a(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P, wl, S0,S1,T0,T1
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        wl = const_pi/2D0 *
     .    (-1D0/3D0*S0**2 + 14D0/15D0*S0*S1 - 9D0/35D0*S1**2
     .     - 0.5D0*T0**2 +1.4D0*T0*T1 + 27D0/70D0*T1**2)
        he3_scatt_l1a = 1D0 + 2D0 * wl / he3_crsect_w(P)
      end

!> Scattering parameter $\lambda_2$ ($\lambda_2^+$)
!> l2 = 1 - 3 &lt;W*sin^4(th/2)sin^2(phi)&gt;/&lt;W&gt;
      function he3_scatt_l2(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P, wl, S0,S1,T0,T1
        call he3_s0s1t0t1(P, S0,S1,T0,T1)
        wl = const_pi/2D0 *
     .    (4D0/15D0*S0**2 - 8D0/21D0*S0*S1 + 52D0/315D0*S1**2
     .     + 0.2D0*T0**2 - 2D0/7D0*T0*T1 + 13D0/105D0*T1**2)
        he3_scatt_l2 = 1D0 - 3D0 * wl / he3_crsect_w(P)
      end

!> Scattering parameter $\gamma_0$,
      function he3_scatt_g0(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_scatt_g0 =
     .    he3_crsect_wi(P) / he3_crsect_w(P)
      end

!> Scattering parameter $\delta_0$,
      function he3_scatt_d0(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_scatt_d0 =
     .    he3_crsect_wd(P) / he3_crsect_w(P)
      end

!><br><img src="img/he3_crsect_w.png">


!> Normal phase viscosity correction c(\lambda), (Sykes-1970)
      function he3_sykes_c(l) !F>
        implicit none
        include 'he3.fh'
        real*8 l,s,ds
        integer n, k
        if (l.eq.1D0) then
          he3_sykes_c = 0.75D0
          return
        endif
        if (l.gt.1D0) then
          he3_sykes_c = NaN
          return
        endif
        s = 0D0
        do n=0,100
          k = (n+1)*(2*n+1)
          ds = (4D0*n + 3D0)/(k*(k-l))
          s = s + ds
          if (ds.lt.1D-8*s) then
            goto 121
          endif
        enddo
121     he3_sykes_c = (1D0-l)/4D0 * s
      end

!> Normal phase viscosity * T^2 (Sykes-1970)
      function he3n_visc(p) !F>
        implicit none
        include 'he3.fh'
        real*8 p,l2,B
        l2 = he3_scatt_l2(p)
        B = 8D0 * const_pi**4 * const_hbar**6
     .    / (he3_meff(p)**3 * const_kb**2 * he3_crsect_w(p))
        he3n_visc = 0.2D0*he3_rho(p)*he3_vf(p)**2*(1D0+he3_f1a/3D0)
     .     * 2D0/const_pi * B/(1D0-l2) * he3_sykes_c(l2)
      end

!> Normal phase thermal conductivity correction H(\lambda), (Sykes-1970)
      function he3_sykes_h(l) !F>
        implicit none
        include 'he3.fh'
        real*8 l,s,ds
        integer n, k
        if (l.eq.3D0) then
          he3_sykes_h = 5D0/12D0
          return
        endif
        if (l.gt.3D0) then
          he3_sykes_h = NaN
          return
        endif
        s = 0D0
        do n=0,100
          k = (n+1)*(2*n+3)
          ds = (4D0*n + 5D0)/(k*(k-l))
          s = s + ds
          if (ds.lt.1D-8*s) then
            goto 122
          endif
        enddo
122     he3_sykes_h = (3D0-l)/4D0 * s
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Quasiparticle lifetimes

!> Normal state quasiparticle lifetime at the Fermi level $\tau_N(0,T)$, s
!> Einzel JLTP32 (1978) p.28,34
!> <br>Einzel JLTP84 (1991) f.4
!> <br>In VW2.38. tau_n0 is different by pi/4 factor!
      function he3_tau_n0(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        real*8 S0,S1,T0,T1, W
        W = he3_crsect_w(p)
        he3_tau_n0 = 32D0 *
     .    he3_tfeff(P)*const_hbar/const_pi**2
     .    / W / const_kb / (1D-3*ttc*he3_tc(P))**2
      end

!> Thermal average of normal state quasiparticle lifetime $3/4\tau_N(0,T)$, s
!> Einzel JLTP84 (1991) f.5
      function he3_tau_n_av(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        he3_tau_n_av = 0.75D0 * he3_tau_n0(ttc, p)
      end

!> Thermal average of normal state spin diffusion transport time, s
!> Einzel JLTP84 (1991) p.328
      function he3_tau_nd(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        he3_tau_nd = 0.75D0 * he3_tau_n0(ttc, p)
     .    / (1D0-he3_scatt_l1a(p))
      end

!> Thermal average ot normal state viscous transport time, s
!> See Einzel JLTP84 (1990) p.41
      function he3_tau_nv(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        he3_tau_nv = 0.75D0 * he3_tau_n0(ttc, p)
     .    / (1D0-he3_scatt_l2(p))
      end

!> Hydrodynamic spin diffusion in normal liquid, cm2/s
!> Einzel JLTP84 (1991) f.23
!> <br>1/3 * vf^2 * tau_n0 * (1+f0a) * 3/4 1/(1-L1)  # Einzel-1991
!> <br>1/3 * vf^2 * tau_n0 * (1+f0a) * f_e(L1)       # VW 2.40 + 2.71
!> <br>Result is same if tau_n0 is different by pi/4 factor
      function he3_diffn_hydr(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, f0a, tau, vf
        tau  = he3_tau_nd(ttc,p)
        vf   = he3_vf(p)
        f0a  = he3_f0a(p)
        he3_diffn_hydr = vf**2 / 3D0 * (1D0+f0a) * tau
      end

!> Frequency-dependent spin diffusion D_perp in normal liquid, cm2/s
!> Einzel JLTP84 (1991) f.22, Bunkov PRL65
      function he3_diffn_perp(ttc, p, nu0) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p, nu0
        real*8 vf, f0a, tau, oe
        f0a  = he3_f0a(p)
        vf   = he3_vf(p)
        tau  = he3_tau_nd(ttc,p)
        oe  = -f0a/(1D0+f0a) * nu0*const_2pi
        he3_diffn_perp = vf**2 / 3D0 * (1D0+f0a)
     .    * tau /(1D0 + (tau * oe)**2)
      end

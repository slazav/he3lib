!HH> He3 Dipole energy and Leggett frequency

!> Dipole coefficient $g_D$, [1/(erg cm^3)] vs P [bar]
!> restored from experimental data by E.Thuneberg.
!> From ROTA texture library.
      function he3_gd(p) !F>
        implicit none
        include 'he3.fh'
        real*8 p
        he3_gd = 1D32 * (0.27733D0 + p*(5.8087D-4 + 2.515D-4*p))
      end

!> Dipole coefficient lambda_d [erg/cm^3] vs T/Tc, P [bar]
!> See Thuneberg-2001 f.5 and f.24.
      function he3_ld(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 p,ttc
        he3_ld = he3_gd(p) * (he3_gap(ttc, p)
     .        * const_kb * 1D-3 * he3_tc(p))**2
      end

!> B-phase Leggett frequency [Hz] vs T/Tc, P[bar]
!> See Thuneberg-2001 f.47.
      function he3_nu_b(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p
        real*8 chi
        chi = he3_chi_b(ttc, p) * he3_chi_n(p)
        he3_nu_b = he3_gyro * dsqrt(15D0 * he3_ld(ttc, p) / chi)
     .   / const_2pi ! rad/s->Hz
      end

!> B-phase Leggett freq, Hz (less accurate formula without use of g_d)
      function he3_nu_b1(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap
        gap  = he3_gap(ttc,p) * const_kb * he3_tc(p)/1D3 ! mk->K
        he3_nu_b1 = dsqrt(3D0 / 8D0 / const_pi /
     .                   he3_chi_b(ttc,p)/he3_chi_n(p))
     .    * he3_gyro**2 * const_hbar * he3_2n0(p) / 2D0
     .    * gap * dlog(he3_tfeff(p)*const_kB/gap)
      end

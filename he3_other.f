!> Extrapolated GL coherence length, cm
!> see Thuneberg-2001, p.667
!> No strong coupling corrections are needed!
      function he3_xigl(ttc,p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,bcsgap
        bcsgap = const_kb*1D-3*he3_tc(p)*he3_bcsgap(ttc)
        he3_xigl = const_hbar * he3_vf(p)
     .    / (dsqrt(10D0)*bcsgap)
      end

!> Equilibrium vortex number
      function he3_vneq(ttc,p,omega,r) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,omega,r
        real*8 kappa,help,rv,rc
        kappa=6.65D-4
        rc=he3_xigl(ttc,p)
        rv=dsqrt(kappa/(const_2pi*omega))
        help=1D0-dsqrt(kappa*LOG(rv/rc)
     .                / (4D0*const_pi*omega*r*r))
!      help=1.0_dp
        he3_vneq = const_2pi*omega*r**2/kappa*help**2
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Leggett-Takagi tau_r [s] vs T/Tc, 20bar
!> Ref: WV pic.10.5 20bar
      function he3_tau_r(ttc) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc
        he3_tau_r = 1.2D-7/dsqrt(1.0D0-ttc)
      end

!> Leggett-Takagi tau_f [s] vs T/Tc, 20bar
      function he3_tau_f(ttc) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, P
        P=20D0
        he3_tau_f = 1D0 /
     .    (const_2pi**2 *he3_nu_b(P,ttc) * he3_tau_r(ttc))
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

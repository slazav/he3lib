! Leggett-Takagi tau_r and tau_s [s] vs T/Tc, 20bar
! Ref: WV pic.10.5 20bar

      function He3_tau_r(TTC)
        implicit none
        include 'he3.fh'
        real*8 TTC
        He3_tau_r = 1.2D-7/dsqrt(1.0D0-TTC)
      end

      function He3_tau_f(TTC)
        implicit none
        include 'he3.fh'
        real*8 TTC, P
        P=20D0
        He3_tau_f = 1D0 /
     .    (4D0*const_pi**2 *He3_Flegg(P, TTC)*He3_tau_r(TTC))
      end

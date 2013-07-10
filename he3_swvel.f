! Osheroff's spin wave velocity. S [cm/s] vs P [bar], T [mK]
! Osheroff's spin wave velocity if recalculated to arbitrary pressure
! by taking value of Osheroff for melting pressure (1100 cm/sek)
! and assuming S is proportional to Fermi velocity.
! Origin: Mukharskii, Dmitriev

      function He3_swvel(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T
        He3_swvel =
     .   1100D0/He3_Vf(34.3D0)*He3_Vf(P)*SQRT(1D0-T/He3_Tc(P))
      end

! Parallel Fomin spin wave velocity. Cpar [cm/c] vs P [bar], T [mK]
! Origin: Mukharskii, Dmitriev

      function he3_swvel_par(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T
        he3_swvel_par=He3_swvel(P,T)*SQRT(2.0D0)
      end

! Perp Fomin spin wave velocity. Cper [cm/c] vs P [bar], T [mK]
! Origin: Mukharskii, Dmitriev

      function he3_swvel_per(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T
        he3_swvel_per=He3_swvel(P,T)*SQRT(1.5D0)
      end


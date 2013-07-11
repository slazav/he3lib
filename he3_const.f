! constants
! Origin: Mukharskii, Dmitriev

      block data he3_const_b
        implicit none
        include 'he3.fh'
        data
     .    he3_Pabn   /21.22D0/,      ! A-B-Normal crit.pt, bar
     .    he3_Tabn   /2.273D0/,      ! A-B-Normal crit.pt, mK
     .    he3_Psmin  /29.3D0/,       ! Melting curve minimum, mK (Daunt, 1959)
     .    he3_Tsmin  /320D0/,        ! Melting curve minimum, mK (Daunt, 1959)
     .    he3_Pa     /34.358D0/,     ! A-N-Solid crit.pt, bar
     .    he3_gyro   /20378.0D0/,    ! g-factor
     .    he3_amass  /5.0079D-24/,   ! He3 atom mass, g
     .    ANA  /6.02214D23/,   ! Avogadro constant, 1/mole
     .    R    /8.314472D+7/,  ! R-gas constant, SGS
     .    HC   /1.05450D-27/,  ! SGS
     .    AKB  /1.38054D-16/,  ! SGS
     .    PI   /3.1415926535897932D0/
      end

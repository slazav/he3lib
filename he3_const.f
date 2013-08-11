      block data he3_const_block
        implicit none
        include 'he3.fh'

        data
     .    he3_gyro   /20378.0D0/,    ! g-factor
     .    he3_amass  /5.0079D-24/,   ! He3 atom mass, g
     .    he3_mmass  /3.016D0/,      ! He3 molar mass, g/mol
     .    const_na   /6.02214129D+23/,  ! Avogadro constant, 1/mole
     .    const_kb   /1.3806488D-16/,   ! SGS
     .    const_r    /8.314472D+7/,     ! R-gas constant, kb*na, SGS
     .    const_h    /6.62606957D-27/,  ! SGS
     .    const_hbar /1.054571726D-27/, ! SGS
     .    const_pi    /3.1415926535897932D0/,
     .    he3_Pcr   /78.111D-3/,
     .    he3_Tcr   /3.324D0/,
     .    he3_Pabn   /21.22D0/,      ! A-B-Normal crit.pt, bar
     .    he3_Tabn   /2.2311D0/,     ! A-B-Normal crit.pt, mK.
                                     ! Graywall -> PLTC correction (2.273->2.2311)
     .    he3_Psmin  /29.3113D0/,    ! Melting curve minimum, bar (PLTS2000)
     .    he3_Tsmin  /0.31524D0/,    ! Melting curve minimum, K (PLTS2000)
     .    he3_Pa     /34.3407D0/,    ! A-N-Solid crit.pt, bar (PLTS2000)
     .    he3_Ta     /2.444D0/,      ! A-N-Solid crit.pt, mK  (PLTS2000)
     .    he3_Pb     /34.3609D0/,    ! A-B-Solid crit.pt, bar (PLTS2000)
     .    he3_Tb     /1.896D0/,      ! A-B-Solid crit.pt, mK  (PLTS2000)
     .    he3_Pneel  /34.3934D0/,    ! A-B-Solid crit.pt, bar (PLTS2000)
     .    he3_Tneel  /0.902D0/       ! A-B-Solid crit.pt, mK  (PLTS2000)

      end
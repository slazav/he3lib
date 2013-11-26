      block data he3_const_block
        implicit none
        include 'he3.fh'

        data
     .    he3_gyro   /20378.0D0/,    ! g-factor
     .    he3_amass  /5.0079D-24/,   ! He3 atom mass, g
     .    he3_mmass  /3.0158281D0/,  ! He3 molar mass, g/mol

     .    const_na   /6.02214129D+23/,  ! Avogadro constant, 1/mole
     .    const_kb   /1.3806488D-16/,   ! SGS
     .    const_r    /8.314472D+7/,     ! R-gas constant, kb*na, SGS
     .    const_h    /6.62606957D-27/,  ! SGS
     .    const_hbar /1.054571726D-27/, ! SGS
     .    const_pi    /3.1415926535897932D0/,
     .    const_2pi   /6.2831853071795864D0/,

     .    he3_Pcr   /1.16317D0/,
     .    he3_Tcr   /3.324D0/,
     .    he3_Pabn   /21.22D0/,      ! A-B-Normal crit.pt, bar
     .    he3_Tabn   /2.273D0/,      ! A-B-Normal crit.pt, mK.

     .    he3_Pm     /29.3113D0/,    ! Melting curve minimum, bar (PLTS2000)
     .    he3_Tm     /0.31524D0/,    ! Melting curve minimum, K (PLTS2000)

     .    he3_Pa     /34.3380D0/,    ! A-N-Solid crit.pt, bar (Greywall-86)
     .    he3_Ta     /2.491D0/,      ! A-N-Solid crit.pt, mK  (Greywall-86)
     .    he3_Pb     /34.3580D0/,    ! A-B-Solid crit.pt, bar (Greywall-86)
     .    he3_Tb     /1.932D0/,      ! A-B-Solid crit.pt, mK  (Greywall-86)
     .    he3_Ps     /34.3905D0/,    ! AFM transition in solid He3, bar (Greywall-86)
     .    he3_Ts     /0.9291D0/,     ! AFM transition in solid He3, mK  (Greywall-86)

     .    he3_Pa_plts  /34.3407D0/,    ! A-N-Solid crit.pt, bar (PLTS2000)
     .    he3_Ta_plts  /2.444D0/,      ! A-N-Solid crit.pt, mK  (PLTS2000)
     .    he3_Pb_plts  /34.3609D0/,    ! A-B-Solid crit.pt, bar (PLTS2000)
     .    he3_Tb_plts  /1.896D0/,      ! A-B-Solid crit.pt, mK  (PLTS2000)
     .    he3_Ps_plts  /34.3934D0/,    ! AFM transition in solid He3, bar (PLTS2000)
     .    he3_Ts_plts  /0.902D0/       ! AFM transition in solid He3, mK  (PLTS2000)

      end
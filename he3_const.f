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

     .    he3_Pm     /29.3113D0/,    ! Melting curve minimum, bar
     .    he3_Tm     /0.31524D0/,    ! Melting curve minimum, K

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
     .    he3_Ts_plts  /0.902D0/,      ! AFM transition in solid He3, mK  (PLTS2000)

     .    rota_rcell    /0.2925D0/,       ! cell radius (2011-2014)
     .    rota_nmra     /96.69139692D0/,  ! field/current in nmrA solenoid [G/A] (normal phase measurements 2014-10-10: f0=832803.65Hz I0=2655.663813mA)
     .    rota_nmrb     /136.6058277D0/,  ! field/current in nmrB solenoid [G/A] (normal phase measurements 2014-10-29, 2014-11-03: f0=1176586.91Hz I0=1907.345mA)
     .    rota_hmina_r  /1.032D0/,        ! effective radius of the HminA coil [cm]
     .    rota_hmina    /2.239D0/,        ! field/current in the center of HminA coil [G/A]
     .    rota_hmina_mr /1.652D0/,        ! quadratic radial term of the HminA field, [G/A/cm^2]
     .    rota_hmina_i0 /9.304D-8/,       ! effectve HminA coil current caused by NMR field distortion, divided by NMR freq [G/Hz]
     .    rota_rrda     /1.0410D-4/       ! radiation damping constant for the nmrA spectrometer

      end
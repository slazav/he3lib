      block data he3_const_block
        implicit none
        include 'he3.fh'

        data
     .    he3_gyro   /20378.0D0/,    ! g-factor, [1/s/G]
     .    he3_amass  /5.0079D-24/,   ! He3 atom mass, [g]
     .    he3_mmass  /3.0158281D0/,  ! He3 molar mass, [g/mol]

     .    const_na   /6.02214129D+23/,  ! Avogadro constant, [1/mole]
     .    const_kb   /1.3806488D-16/,   ! Boltzmann constant, [erg/K]
     .    const_r    /8.314472D+7/,     ! R-gas constant, kb*na, [erg/K/mol]
     .    const_h    /6.62606957D-27/,  ! Planck constant, [g*cm2/s]
     .    const_hbar /1.054571726D-27/, ! reduced Planck constant, [g*cm2/s]
     .    const_mu0  /1.2566370614D0/,  ! vacuum permeability [G*cm/A]
     .    const_ev   /1.602176634D-12/, ! electronvolt [erg]
     .    const_e     /2.7182818284590452D0/, ! e
     .    const_pi    /3.1415926535897932D0/, ! pi
     .    const_2pi   /6.2831853071795864D0/, ! 2*pi
     .    const_euler /0.5772156649015329D0/, ! Euler's constant
     .    const_z2    /1.6449340668482264D0/, ! zeta(2)
     .    const_z3    /1.2020569031595943D0/, ! zeta(3)
     .    const_z4    /1.0823232337111382D0/, ! zeta(4)
     .    const_z5    /1.0369277551433699D0/, ! zeta(5)

     .    he3_Pcr   /1.16317D0/,
     .    he3_Tcr   /3.324D0/,
     .    he3_Pabn   /21.22D0/,      ! A-B-Normal crit.pt, bar
     .    he3_Tabn   /2.273D0/,      ! A-B-Normal crit.pt, mK.

     .    he3_Pm     /29.3113D0/,    ! Melting curve minimum, bar (PLTS2000)
     .    he3_Tm     /0.31524D0/,    ! Melting curve minimum, K   (PLTS2000)

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

     .    he3_Pabn_plts /21.222D0/,    ! A-B-Normal pt., Greywall->PLTS, bar
     .    he3_Tabn_plts /2.2315D0/,    ! A-B-Normal pt., Greywall->PLTS, mK

     .    rota_rcell    /0.2925D0/,       ! cell radius (2011-2014)
     .    rota_nmra     /96.69139692D0/,  ! field/current in nmrA solenoid [G/A] (normal phase measurements 2014-10-10: f0=832803.65Hz I0=2655.663813mA)
     .    rota_nmrb     /136.6058277D0/,  ! field/current in nmrB solenoid [G/A] (normal phase measurements 2014-10-29, 2014-11-03: f0=1176586.91Hz I0=1907.345mA)
     .    rota_hmina_r  /1.032D0/,        ! effective radius of the HminA coil [cm]
     .    rota_hmina_n  /4D0/,            ! number of turns of the HminA coil
     .    rota_hmina    /2.2305D0/,       ! field/current in the center of HminA coil [G/A]
     .    rota_hmina_mr /1.652D0/,        ! quadratic radial term of the HminA field, [G/A/cm^2]
     .    rota_hmina_i0i /-0.02918D0/,    ! effectve HminA coil current divided by NMR current
     .    rota_hmina_i0f /-9.3050D-08/,   ! effectve HminA coil current divided by NMR frequency
                                          ! rota_hmina_i0i = rota_hmina_i0f * he3_gyro/2/pi * rota_nmra
     .    rota_rrda     /1.01D-4/         ! radiation damping constant for the nmrA spectrometer

      end

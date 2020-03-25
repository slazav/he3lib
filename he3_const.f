!H> Constants

      block data he3_const_block
        implicit none
        include 'he3.fh'

        data
!H> Mathematical constants
     .    const_e     /2.7182818284590452D0/, !C> e
     .    const_pi    /3.1415926535897932D0/, !C> pi
     .    const_2pi   /6.2831853071795864D0/, !C> 2*pi
     .    const_euler /0.5772156649015329D0/, !C> Euler's constant
     .    const_z2    /1.6449340668482264D0/, !C> zeta(2)
     .    const_z3    /1.2020569031595943D0/, !C> zeta(3)
     .    const_z4    /1.0823232337111382D0/, !C> zeta(4)
     .    const_z5    /1.0369277551433699D0/, !C> zeta(5)

!H> Physical constants
     .    const_na   /6.02214129D+23/,  !C> Avogadro constant, [1/mole]
     .    const_kb   /1.3806488D-16/,   !C> Boltzmann constant, [erg/K]
     .    const_r    /8.314472D+7/,     !C> R-gas constant, kb*na, [erg/K/mol]
     .    const_h    /6.62606957D-27/,  !C> Planck constant, [g*cm2/s]
     .    const_hbar /1.054571726D-27/, !C> reduced Planck constant, [g*cm2/s]
     .    const_mu0  /1.2566370614D0/,  !C> vacuum permeability [G*cm/A]
     .    const_ev   /1.602176634D-12/, !C> electronvolt [erg]

     .    he3_gyro   /20378.0D0/,    !C> He3: g-factor, [1/s/G]
     .    he3_amass  /5.0079D-24/,   !C> He3: atom mass, [g]
     .    he3_mmass  /3.0158281D0/,  !C> He3: molar mass, [g/mol]

!H> Helium-3 phase diagram constants
     .    he3_pcr   /1.16317D0/,     !C> He3: critical pressure [bar]
     .    he3_tcr   /3.324D0/,       !C> He3: critical temperature [K]
     .    he3_pabn   /21.22D0/,      !C> He3: A-B-Normal tricritical pressure [bar]
     .    he3_tabn   /2.273D0/,      !C> He3: A-B-Normal tricritical temperature [mK]

     .    he3_pm     /29.3113D0/,    !C> He3: Melting curve minimum [bar] (PLTS2000)
     .    he3_tm     /0.31524D0/,    !C> He3: Melting curve minimum [K]   (PLTS2000)

     .    he3_pa     /34.3380D0/,    !C> He3: A-N-Solid crit.pt, bar (Greywall-86)
     .    he3_ta     /2.491D0/,      !C> He3: A-N-Solid crit.pt, mK  (Greywall-86)
     .    he3_pb     /34.3580D0/,    !C> He3: A-B-Solid crit.pt, bar (Greywall-86)
     .    he3_tb     /1.932D0/,      !C> He3: A-B-Solid crit.pt, mK  (Greywall-86)
     .    he3_ps     /34.3905D0/,    !C> Solid He3: AFM transition at melting curve [bar] (Greywall-86)
     .    he3_ts     /0.9291D0/,     !C> Solid He3: AFM transition at melting curve [mK]  (Greywall-86)

!H> Helium-3 phase diagram constants (PLTS2000)
     .    he3_pa_plts  /34.3407D0/,  !C> He3: A-N-Solid crit.pt, bar (PLTS2000)
     .    he3_ta_plts  /2.444D0/,    !C> He3: A-N-Solid crit.pt, mK  (PLTS2000)
     .    he3_pb_plts  /34.3609D0/,  !C> He3: A-B-Solid crit.pt, bar (PLTS2000)
     .    he3_tb_plts  /1.896D0/,    !C> He3: A-B-Solid crit.pt, mK  (PLTS2000)
     .    he3_ps_plts  /34.3934D0/,  !C> He3: AFM transition in solid He3, bar (PLTS2000)
     .    he3_ts_plts  /0.902D0/,    !C> He3: AFM transition in solid He3, mK  (PLTS2000)

     .    he3_pabn_plts /21.222D0/,  !C> He3: A-B-Normal pt., Greywall->PLTS, bar
     .    he3_tabn_plts /2.2315D0/,  !C> He3: A-B-Normal pt., Greywall->PLTS, mK

!H> ROTA-specific constants
     .    rota_rcell    /0.2925D0/,       !C> ROTA: cell radius (2011-2014)
     .    rota_nmra     /96.69139692D0/,  !C> ROTA: field/current in nmrA solenoid [G/A] (normal phase measurements 2014-10-10: f0=832803.65Hz I0=2655.663813mA)
     .    rota_nmrb     /136.6058277D0/,  !C> ROTA: field/current in nmrB solenoid [G/A] (normal phase measurements 2014-10-29, 2014-11-03: f0=1176586.91Hz I0=1907.345mA)
     .    rota_hmina_r  /1.032D0/,        !C> ROTA: effective radius of the HminA coil [cm]
     .    rota_hmina_n  /4D0/,            !C> ROTA: number of turns of the HminA coil
     .    rota_hmina    /2.2305D0/,       !C> ROTA: field/current in the center of HminA coil [G/A]
     .    rota_hmina_mr /1.652D0/,        !C> ROTA: quadratic radial term of the HminA field, [G/A/cm^2]
     .    rota_hmina_i0i /-0.02918D0/,    !C> ROTA: effectve HminA coil current divided by NMR current
     .    rota_hmina_i0f /-9.3050D-08/,   !C> ROTA: effectve HminA coil current divided by NMR frequency, rota_hmina_i0i = rota_hmina_i0f * he3_gyro/2/pi * rota_nmra
     .    rota_rrda     /1.01D-4/         !C> ROTA: radiation damping constant for the nmrA spectrometer

      end

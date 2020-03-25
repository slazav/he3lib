!HH> Constants

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

     .    h1_gyro    /26752.218744D0/,  !C> H1 g-factor, [rad/s/G]
     .    h2_gyro    /4106.5D0/,        !C> H2 g-factor, [rad/s/G]
     .    h_amass    /1.67355755D-24/,  !C> H molar mass, [g/mol]
     .    h_mmass    /1.00784D0/        !C> H molar mass, [g/mol]

     .    he3_gyro   /20378.9D0/,       !C> He3 g-factor, [rad/s/G]
     .    he3_amass  /5.00789994D-24/,  !C> He3 atom mass, [g]
     .    he3_mmass  /3.0158281D0/,     !C> He3 molar mass, [g/mol]

     .    he4_amass  /6.64647641D-24/,  !C> He4 atom mass, [g]
     .    he4_mmass  /4.002602D0/       !C> He4 molar mass, [g/mol]

      end

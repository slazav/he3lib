!H> Helium-3 as a normal fermi-liquid

!> See Wheatley, Rev.Mod.Phys. 47, 415(1975) /tables at p.467/

!H> Molar volume and specific heat

!>    He3 molar volume [cm^3/mole] (exp data, Greywall-86)
      function he3_vm(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
!          he3_vm = ! Graywall-83
!     .     +6.218603D-08 * P**6
!     .     -7.287781D-06 * P**5
!     .     +3.393465D-04 * P**4
!     .     -8.190787D-03 * P**3
!     .     +1.171477D-01 * P**2
!     .     -1.269684D+00 * P**1
!     .     +3.683605D+01
          he3_vm =   ! Graywall-86 (from Wheatley-75)
     .       + 36.837231D0
     .       - 0.11803474D+1*P
     .       + 0.83421417D-1*P**2
     .       - 0.38859562D-2*P**3
     .       + 0.94759780D-4*P**4
     .       - 0.91253577D-6*P**5
        else
          he3_vm = NaN
        endif
      end

!>    He3 specific heat Cv/RT  [1/K], (exp data, Greywall-86)
!>    see also Alvesalo PRL44 1076 (1980) - they have different values!
      function he3_gammaf(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
          he3_gammaf =
     .     +0.27840464D+1
     .     +0.69575243D-1*P
     .     -0.14738303D-2*P**2
     .     +0.46153498D-4*P**3
     .     -0.53785385D-6*P**4
        else
          he3_gammaf = NaN
        endif
      end

!H> derived values

!> He3 heat capacity (C/R) vs T(K) and P(bar)
      function he3_c_n(T,P) !F>
        implicit none
        include 'he3.fh'
        real*8 T,P
        he3_c_n = he3_gammaf(P)*T
      end

!> He3 density [g/cm^3]
      function he3_rho(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_rho = he3_mmass / he3_vm(P)
      end

!> 2N0 [1/erg/cm^3]
      function he3_2n0(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
!        he3_2n0 = he3_gammaf(P) / he3_vm(P) / 7.5421D-40
        he3_2n0 = he3_gammaf(P) / he3_vm(P) *
     .    3D0 *const_na/const_kb/const_pi**2
      end

!> He3 Fermi momentum [g cm/s] vs P [bar]
      function he3_pf(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
!        he3_pf = 2.7551D-19 / he3_vm(P)**.3333333D0
        he3_pf = const_h * (3D0/8D0/const_pi *
     .    const_na/he3_vm(P))**.3333333D0
      end

!> He3 Fermi velocity [cm/s] vs P [bar]
      function he3_vf(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
!        he3_vf = 4946.6423D0 * he3_vm(P)**.3333333D0 / he3_gammaf(P)
        he3_vf = he3_pf(P)/he3_meff(P)
      end

!> He3 effective mass m_eff/m_3 vs P [bar]
      function he3_mm(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
!        he3_mm = 11.1217D0 * he3_gammaf(P) / he3_vm(P)**0.6666667D0
        he3_mm = const_h**3/8D0/const_pi *
     .           he3_2n0(P)/he3_pf(P) / he3_amass
      end

!> He3 effective mass [g] vs P [bar].
      function he3_meff(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
!        he3_meff = 5.569812141D-23 * he3_gammaf(P) / he3_vm(P)**0.6666667D0
        he3_meff = const_h**3/8D0/const_pi *
     .             he3_2n0(P)/he3_pf(P)
      end

!> He3 F1s fermi-liquid parameter
      function he3_f1s(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_f1s = 3D0*(he3_mm(P)-1D0)
      end

!> He3 average atomic spacing [angstr]
      function he3_a(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_a = (he3_vm(P)/const_na)**0.3333333D0 * 1D8
      end

!> He3 average dipolar coupling energy, K
      function he3_gdk(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_gdk = 2D0/3D0*const_pi*he3_gyro**2 / he3_vm(P)
     .          * const_hbar**2 * const_na/const_kb
      end

!> He3 effective Fermi temperature, K
      function he3_tfeff(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_tfeff = const_pi**2/2D0 / he3_gammaf(P)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H>  Sound velocity and F0s parameter

!> First sound velocity c1 (exp data), m/s, from Wheatley-75
      function he3_c1(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
          he3_c1 =
     .     -4.604822D-07*P**6
     .     +5.472623D-05*P**5
     .     -2.635920D-03*P**4
     .     +6.846620D-02*P**3
     .     -1.130160D+00*P**2
     .     +1.764852D+01*P
     .     +1.829327D+02
          he3_c1 = he3_c1 * 1D2 ! m/s -> cm/s
        else
          he3_c1 = NaN
        endif
      end

!> He3 F0s fermi-liquid parameter
      function he3_f0s(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_f0s = 3D0*he3_amass*he3_meff(P)
     .             *he3_c1(P)**2 / he3_pf(P)**2 - 1D0
      end

!>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Susceptibility and F0a parameter

!> He3 F0a fermi-liquid parameter (same as Z0/4)
!> Wheatley-75
!> Ramm, JLTP 2 539 (1970)
!> + Hensley, JLTP89 501 (1992), JLTP90 149 (1993)
      function he3_f0a(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P, tmag
        if (P.ge.0D0.and.P.le.34.40D0) then
!!!        Magnetic temperature T* (exp data), K, from Wheatley-75
!          tmag =
!     .     +2.119835D-09*P**6
!     .     -2.382702D-07*P**5
!     .     +1.043133D-05*P**4
!     .     -2.283417D-04*P**3
!     .     +2.794162D-03*P**2
!     .     -2.427311D-02*P
!     .     +3.588074D-01
!        he3_f0a = 
!     .   3D0*const_kb*tmag(P)*he3_meff(P)/he3_pf(P)**2 - 1D0
!!!       Hensley JLTP90 149 (1993)
        he3_f0a =
     .    +1.0891D-09*P**6
     .    -1.3033D-07*P**5
     .    +6.0659D-06*P**4
     .    -1.3904D-04*P**3
     .    +1.6950D-03*P**2
     .    -1.2308D-02*P   
     .    -6.9863D-01
        else
          he3_f0a = NaN
        endif
      end

!> Susceptibility [sgs] vs P [bar], T [mK]
!> see Einzel-1991 f.10
      function he3_chi_n(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        he3_chi_n = he3_2n0(P)*(he3_gyro*const_hbar/2D0)**2
     .    / (1D0 + he3_f0a(P))
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Other fermi-liquid parameters

!> He3 F1a fermi-liquid parameter (Zavjalov-2015, from spin-wave velocities)
!> Corruccini PRL27 650 (1971) -- Leggett-Rice effect in 3He-N, not accurate
!> Osheroff PhB90 20 (1977) -- Spin-wave velocity in 3He-B, not accurate
!> Greywall-1983 -- high temperature Cv
!> Zavjalov-2015 -- Spin-wave velocity in 3He-B
!> theory, spin waves: Dorfle PRB23 3267 (1981) + F3s
!> theory, spin waves: Cross JLTP 21 525 (1975)
      function he3_f1a(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
! ! Greywall-83
!          he3_f1a =
!     .      +1.489333D-08 * P**6
!     .      -1.661179D-06 * P**5
!     .      +7.056909D-05 * P**4
!     .      -1.430616D-03 * P**3
!     .      +1.476252D-02 * P**2
!     .      -9.170753D-02 * P**1
!     .      -5.506076D-01
!  our data
          he3_f1a = -0.598D0 -0.00214D0*P
        else
          he3_f1a = NaN
        endif
      end

!> He3 F2a fermi-liquid parameter (zero at the moment)
!> Halperin???
      function he3_f2a(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
          he3_f2a = 0D0
        else
          he3_f2a = NaN
        endif
      end

!> He3 F2s fermi-liquid parameter (zero at the moment)
!>  Engel, Ihas, Phys. Rev. Lett. 55, 955958 (1985)
!>  Hamot, Lee, ... Halperin, JLTP 99 p651 (1995)
!>  Mastumoto et al. JLTP 102 p227 (1996)
      function he3_f2s(P) !F>
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
          he3_f2s = 0D0
!          ! value from Engel-1985
!          he3_f2s = -1.264D0 + 0.896D0*dsqrt(P)
!     .              -0.187D0*P + 0.0163D0*dsqrt(P**3)
!          ! value from Mastumoto-1996
!          he3_f2s = 0.15D0;
        else
          he3_f2s = NaN
        endif
      end



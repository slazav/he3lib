! Fermi-liquid parameters
! See Wheatley, Rev.Mod.Phys. 47, 415(1975) /tables at p.467/

!     Molar volume (exp data)
      function He3_Vm(P) ! Greywall-83
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
!          He3_Vm = ! Graywall-83
!     .     +6.218603D-08 * P**6
!     .     -7.287781D-06 * P**5
!     .     +3.393465D-04 * P**4
!     .     -8.190787D-03 * P**3
!     .     +1.171477D-01 * P**2
!     .     -1.269684D+00 * P**1
!     .     +3.683605D+01
          He3_Vm =   ! Graywall-86 (from Wheatley-75)
     .       + 36.837231D0
     .       - 0.11803474D+1*P
     .       + 0.83421417D-1*P**2
     .       - 0.38859562D-2*P**3
     .       + 0.94759780D-4*P**4
     .       - 0.91253577D-6*P**5
        else
          He3_Vm = NaN
        endif
      end

!     C_p/RT (exp data), 1/(mol K), Greywall-86
      function He3_gammaf(P)
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
          He3_gammaf =
     .     +0.27840464D+1
     .     +0.69575243D-1*P
     .     -0.14738303D-2*P**2
     .     +0.46153498D-4*P**3
     .     -0.53785385D-6*P**4
        else
          He3_gammaf = NaN
        endif
      end

!     First sound velosity c1 (exp data), m/s, from Wheatley-75
      function He3_c1(P)
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
          He3_c1 =
     .     -4.604822D-07*P**6
     .     +5.472623D-05*P**5
     .     -2.635920D-03*P**4
     .     +6.846620D-02*P**3
     .     -1.130160D+00*P**2
     .     +1.764852D+01*P
     .     +1.829327D+02
          He3_c1 = He3_c1 * 1D2 ! m/s -> cm/s
        else
          He3_c1 = NaN
        endif
      end

!     Magnetic temperature T* (exp data), K, from Wheatley-75
      function He3_tmag(P)
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
          He3_tmag =
     .     +2.119835D-09*P**6
     .     -2.382702D-07*P**5
     .     +1.043133D-05*P**4
     .     -2.283417D-04*P**3
     .     +2.794162D-03*P**2
     .     -2.427311D-02*P
     .     +3.588074D-01
        else
          He3_tmag = NaN
        endif
      end

!!!!!!!
!  Engel, Ihas, Phys. Rev. Lett. 55, 955958 (1985)
!  Hamot, Lee, ... Halperin, JLTP 99 p651 (1995)



!!!!!!!    derived values

!     density, g/cm^3
      function He3_rho(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_rho = he3_mmass / he3_vm(P)
      end

!     2N0
      function He3_2N0(P)
        implicit none
        include 'he3.fh'
        real*8 P
!        He3_2N0 = he3_gammaf(P) / he3_vm(P) / 7.5421D-40
        He3_2N0 = he3_gammaf(P) / he3_vm(P) *
     .    3D0 *const_na/const_kb/const_pi**2
      end

! Fermi momentum [sgs] vs P [bar]
      function He3_Pf(P)
        implicit none
        include 'he3.fh'
        real*8 P
!        He3_Pf = 2.7551D-19 / He3_Vm(P)**.3333333D0
        He3_Pf = const_h * (3D0/8D0/const_pi *
     .    const_na/He3_Vm(P))**.3333333D0
      end

! Fermi velocity [cm/s] vs P [bar]
      function He3_Vf(P)
        implicit none
        include 'he3.fh'
        real*8 P
!        He3_Vf = 4946.6423D0 * He3_Vm(P)**.3333333D0 / he3_gammaf(P)
        He3_Vf = he3_pf(P)/he3_meff(P)
      end

! Effective mass m_eff/m_3 vs P [bar]
      function He3_mm(P)
        implicit none
        include 'he3.fh'	
        real*8 P
!        He3_mm = 11.1217D0 * he3_gammaf(P) / he3_vm(P)**0.6666667D0
        He3_mm = const_h**3/8D0/const_pi *
     .           he3_2n0(P)/he3_pf(P) / he3_amass
      end

! Effective mass [g] vs P [bar].
      function He3_Meff(P)
        implicit none
        include 'he3.fh'
        real*8 P
!        He3_Meff = 5.569812141D-23 * he3_gammaf(P) / he3_vm(P)**0.6666667D0
        He3_Meff = const_h**3/8D0/const_pi *
     .             he3_2n0(P)/he3_pf(P)
      end

! Suseptibility [sgs] vs P [bar], T [mK]
! Origin: Mukharskii, Dmitriev program
      function He3_chi_n(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_chi_n = 0.25D0 * he3_gyro**2 * const_hbar**2
     .   * const_na/const_kb**2/const_pi**2 * He3_gammaf(P)
     .   / he3_tmag(P) / he3_meff(P) * he3_pf(P)**2 / He3_Vm(P)
      end

      function He3_F0s(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_F0s = 3D0*he3_amass*he3_meff(P)
     .             *he3_c1(P)**2 / he3_pf(P)**2 - 1D0
      end

      function He3_F1s(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_F1s = 3D0*(he3_mm(P)-1D0)
      end

      ! see also Halperin, JLTP 89 (1982)
      ! same as Z0/4
      function He3_F0a(P) ! Greywall-83
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
        he3_f0a = 
     .   3D0*const_kb*he3_tmag(P)*he3_meff(P)/he3_pf(P)**2 - 1D0
        else
          He3_F0a = NaN
        endif
      end

      function He3_F1a(P) ! Greywall-83
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.34.40D0) then
          He3_F1a =
     .      +1.489333D-08 * P**6
     .      -1.661179D-06 * P**5
     .      +7.056909D-05 * P**4
     .      -1.430616D-03 * P**3
     .      +1.476252D-02 * P**2
     .      -9.170753D-02 * P**1
     .      -5.506076D-01
        else
          He3_F1a = NaN
        endif
      end

!     average atomic spacing, angstr.
      function He3_a(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_a = (he3_vm(P)/const_na)**0.3333333D0 * 1D8
      end

!     average dipolar coupling energy, K
      function He3_gdk(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_gdk = 2D0/3D0*const_pi*he3_gyro**2 / he3_vm(P)
     .          * const_hbar**2 * const_na/const_kb
      end

!     effective fermi temperature, K
      function He3_tfeff(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_tfeff = const_pi**2/2D0 / he3_gammaf(P)
      end

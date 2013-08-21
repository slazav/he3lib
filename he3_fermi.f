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
! Einzel-1991 f.10
      function He3_chi_n(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_chi_n = he3_2n0(P)/4D0 *(he3_gyro*const_hbar/2D0)**2
     .    / (1D0 + he3_f0a(P))
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

!     heat capacity, Cv
!     Greywall-1983
      function He3_cv_n(t, v)
        implicit none
        include 'he3.fh'
        real*8 t,v,a,b,c,d
        real*8 s1,s2,s3
        integer i,j
        dimension a(5,4), b(4,3), c(3,3), d(3)
        a(1,1) = -2.9190414D0
        a(1,2) =  5.2893401D+2
        a(1,3) = -1.8869641D+4
        a(1,4) =  2.6031315D+5
        a(2,1) =  0D0
        a(2,2) =  0D0
        a(2,3) =  0D0
        a(2,4) =  0D0
        a(3,1) = -2.4752597D+3
        a(3,2) =  1.8377260D+5
        a(3,3) = -3.4946553D+6
        a(3,4) =  0D0
        a(4,1) =  3.8887481D+4
        a(4,2) = -2.8649769D+6
        a(4,3) =  5.2526785D+7
        a(4,4) =  0D0
        a(5,1) = -1.7505655D+5
        a(5,2) =  1.2809001D+7
        a(5,3) = -2.3037701D+8
        a(5,4) =  0D0

        b(1,1) = -6.5521193D-2
        b(1,2) =  1.3502371D-2
        b(1,3) =  0D0
        b(2,1) =  4.1359033D-2
        b(2,2) =  3.8233755D-4
        b(2,3) = -5.3468396D-5
        b(3,1) =  5.7976786D-3
        b(3,2) = -6.5611532D-4
        b(3,3) =  1.2689707D-5
        b(4,1) = -3.8374623D-4
        b(4,2) =  3.2072581D-5
        b(4,3) = -5.3038906D-7

        c(1,1) = -2.5482958D+1
        c(1,2) =  1.6416936D+0
        c(1,3) = -1.5110378D-2
        c(2,1) =  3.7882751D+1
        c(2,2) = -2.8769188D+0
        c(2,3) =  3.5751181D-2
        c(3,1) =  2.4412956D+1
        c(3,2) = -2.4244083D+0
        c(3,3) =  6.7775905D-2

        d(1) = -7.1613436D+0
        d(2) =  6.0525139D-1
        d(3) = -7.1295855D-3

        if (t.lt.0.1) then
          s1=0D0
          do i=1,5
            do j=0,3
              s1 = s1 + a(i,j+1) * v**(-j) * t**i
            enddo
          enddo
          he3_cv_n = s1
          return
        endif

        if (t.ge.0.1.and.t.lt.2.5) then
          s1=0D0
          do i=0,3
            do j=0,2
              s1 = s1 + b(i+1,j+1) * v**j * t**(-i)
            enddo
          enddo
          s2=0D0
          do i=1,3
            do j=0,2
              s2 = s2 + c(i,j+1) * v**j * t**(-i)
            enddo
          enddo
          s3=0D0
          do j=0,2
            s3 = s3 + d(j+1) * v**j
          enddo
          he3_cv_n = s1 + dexp(-s3/t) * s2
          return
        endif
        he3_cv_n = NaN
      end




!!! He3 phase diagram

! Vapor pressure [bars] vs T [K]
! Sherman, R.H.; Sydoriak, S.G.; Roberts, T.R.
! The 1962 He3 scale of temperatures
! http://archive.org/details/jresv68An6p579

      function He3_Pvap(T)
        implicit none
        include 'he3.fh'
        real*8 T
        if (T.gt.0.2D0.and.T.le.he3_tcr) then
          He3_Pvap =
     .     - 2.49174D0 / T
     .     + 4.80386D0
     .     - 0.286001D0 * T
     .     + 0.198608D0 * T**2
     .     - 0.0502237D0 * T**3
     .     + 0.00505486D0 * T**4
          He3_Pvap = dexp(He3_Pvap) * 1.333224D-3 ! torr -> bar
        else
          He3_Pvap = NaN
        endif
        return
      end

! Melting pressure [bars] vs T [mK]
! Arg: T = 0.0009 .. 31 [K]
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520 (1985)
! Ref: Osborne, Abraham, Weinstock, 1951, 1952
! Ref: Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955)
! see also: Johnson, Symko, Weatley -- Phis. Rev. Lett. 23, 1017 (1969)
      function He3_Pmelt(T)
        implicit none
        include 'he3.fh'
        real*8 T
        if (T.gt.9D-4.and.T.le.0.25D0) then
!!        PLTS2000
          He3_Pmelt =
     .     - 1.3855442D-12 * T**(-3)
     .     + 4.5557026D-9  * T**(-2)
     .     - 6.4430869D-6  * T**(-1)
     .     + 3.4467434D0
     .     - 4.4176438D0 * T**1
     .     + 1.5417437D1 * T**2
     .     - 3.5789858D1 * T**3
     .     + 7.1499125D1 * T**4
     .     - 1.0414379D2 * T**5
     .     + 1.0518538D2 * T**6
     .     - 6.9443767D1 * T**7
     .     + 2.6833087D1 * T**8
     .     - 4.5875709D0 * T**9
         He3_Pmelt = He3_Pmelt * 10D0 ! MPa -> bar
        else if (T.gt.0.25D0.and.T.le.0.5D0) then
!         Interpolation
          He3_Pmelt = 33.2505D0
     .     - 24.6026D0 *T
     .     + 34.5713D0 *T**2
     .     + 20.4086D0 *T**3
     .     - 26.3258D0 *T**4
        else if (T.gt.0.5D0.and.T.le.1.5D0) then
!         Osborne, Abraham, Weinstock, 1951
!         Pm = 26.8 + 13.1 T^2 [atm], T = 0.5 .. 1.5
!         atm->bar: 1.01325
          He3_Pmelt = (26.8D0 + 13.1D0 * T**2) * 1.01325D0
        else if (T.gt.1.5D0.and.T.le.2.06D0) then
!         Interpolation
          He3_Pmelt = 24.787D0
     .    + 9.4456D0 *T
     .    + 0.0933D0 *T**2
     .    + 7.6155D0 *T**3
     .    - 1.5327D0 *T**4
        else if (T.gt.2.06D0.and.T.le.31D0) then
!         Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955))
!         Pm = 25.16 + 20.08201 T^1.517083  [kg/cm2] P=76-3500
!         He4: -17.80 + 17.31457 T^1.555414 [kg/cm2] P=37-3500
!         kgf/cm2 -> bar: 0.980665
          He3_Pmelt = (25.16D0 + 20.08201D0 * T**1.517083D0)
     .                * 0.980665D0
        else
          He3_Pmelt = NaN
        endif
        return
      end

! Greywall-86 melting pressure temperature scale [bars] vs T [mK]
! Arg: T = 0.0009 .. 31K [K]
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520 (1986)
      function He3_Pmelt_gr(T)
        implicit none
        include 'he3.fh'
        real*8 T
        if (T.gt.9D-4.and.T.le.0.25D0) then
          He3_Pmelt_gr = He3_Pa
     .     - 0.19652970D-1 * (T*1D3)**(-3)
     .     + 0.61880268D-1 * (T*1D3)**(-2)
     .     - 0.78803055D-1 * (T*1D3)**(-1)
     .     + 0.13050600D0
     .     - 0.43519381D-1 * (T*1D3)
     .     + 0.13752791D-3 * (T*1D3)**2
     .     - 0.17180436D-6 * (T*1D3)**3
     .     - 0.22093906D-9 * (T*1D3)**4
     .     + 0.85450245D-12* (T*1D3)**5
        else
          He3_Pmelt_gr = He3_Pmelt(T)
        endif
        return
      end

! T_c [mK] vs P [bar]
! Arg: P = 0 .. Pa [bar]
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520 (1985)
      function He3_Tc(P)
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.ge.0D0.and.P.le.he3_pa) then
          He3_Tc= .92938375D0
     .       + .13867188D0*P
     .       - .69302185D-2*P**2
     .       + .25685169D-3*P**3
     .       - .57248644D-5*P**4
     .       + .53010918D-7*P**5
             He3_Tc = 0.96756D0*He3_Tc + 0.031803D0  ! Greywall -> PLTC temp scale
        else
          He3_Tc = NaN
        endif
        return
      end

! T_ab [mK] vs P [bar]
! Arg: P = 0 .. Pmelt [bar]
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520 (1985)
! Note: Pabn=21.22 bar, Tabn = 2.273 mK
      function He3_Tab(P)
        implicit none
        include 'he3.fh'
        real*8 P, Pr
        Pr = P - 21.22D0
        if (Pr.lt.0D0) then
          He3_Tab=He3_Tc(P)
        else
          He3_Tab= 2.273D0
     .       - .10322623D-1*Pr
     .       - .53633181D-2*Pr**2
     .       + .83437032D-3*Pr**3
     .       - .61709783D-4*Pr**4
     .       + .17038992D-5*Pr**5
          He3_Tab = 0.96756D0*He3_Tab + 0.031803D0  ! Greywall -> PLTC temp scale
        endif
        if (P.lt.0D0.or.P.gt.He3_Pb) then
          He3_Tab = NaN
        endif
        return
      end

! T_ab [mk] vs P, H
! Inseob Hahn PhD thesis, p79
! real*8 p0,p1,p2,q1,pp
! real*8 Bc,f3,f4,Pa
! Pa = 34.338
! pp = P/Pa
! Bc=(3391D0 + 21500D0*pp - 8490D0*pp**2) / (1+2.098D0*pp)
! f3 = 1.41D0
! f4 = -0.29D0 -0.41D0*pp
! f5 = (1-f3*Ta**6 - f4*Ta**8 - (1-f3-f4)*Ta**2 + (1+2f3+3f4)*(Ta**4-Ta**2))/
!        /

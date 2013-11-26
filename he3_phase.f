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
     .     + 2.24846D0 * dlog(T)
          He3_Pvap = dexp(He3_Pvap) * 1.333224D-3 ! torr -> bar
        else
          He3_Pvap = NaN
        endif
        return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Melting pressure [bars] vs T [mK]
! Arg: T = 0.0009 .. 31 [K]
! Ref: Greywall. PRB33 7520 (1986) f.A1
! Ref: Osborne, Abraham, Weinstock, 1951, 1952
! Ref: Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955)
! see also: Johnson, Symko, Weatley -- Phis. Rev. Lett. 23, 1017 (1969)
      function He3_Pmelt(T)
        implicit none
        include 'he3.fh'
        real*8 T
        if (T.gt.9D-4.and.T.le.0.25D0) then
          He3_Pmelt = 34.3380D0
     .     - 0.19652970D-1 * (T*1D3)**(-3)
     .     + 0.61880268D-1 * (T*1D3)**(-2)
     .     - 0.78803055D-1 * (T*1D3)**(-1)
     .     + 0.13050600D0
     .     - 0.43519381D-1 * (T*1D3)
     .     + 0.13752791D-3 * (T*1D3)**2
     .     - 0.17180436D-6 * (T*1D3)**3
     .     - 0.22093906D-9 * (T*1D3)**4
     .     + 0.85450245D-12* (T*1D3)**5
        else if (T.gt.0.25D0.and.T.le.0.5D0) then
!         Interpolation (see tests/pmelt_interp.m)
          He3_Pmelt =
     .      3786.422495D0 * T**5
     .     -6984.120614D0 * T**4
     .     +5051.847883D0 * T**3
     .     -1756.109750D0 * T**2
     .      +289.300690D0 * T
     .       +11.551436D0
        else if (T.gt.0.5D0.and.T.le.1.5D0) then
!         Osborne, Abraham, Weinstock, 1951
!         Pm = 26.8 + 13.1 T^2 [atm], T = 0.5 .. 1.5
!         atm->bar: 1.01325
          He3_Pmelt = (26.8D0 + 13.1D0 * T**2) * 1.01325D0
        else if (T.gt.1.5D0.and.T.le.2D0) then
!         Interpolation (see tests/pmelt_interp.m)
          He3_Pmelt =
     .      -53.992350D0 * T**3
     .     +286.394879D0 * T**2
     .     -454.915548D0 * T
     .     +277.229670D0
        else if (T.gt.2D0.and.T.le.31D0) then
!         Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955))
!         Pm = 25.16 + 20.08201 T^1.517083  [kg/cm2] P=76-3500
!         He4: -17.80 + 17.31457 T^1.555414 [kg/cm2] P=37-3500
!         kgf/cm2 -> bar: 0.980665
          He3_Pmelt = (25.16D0 + 20.08201D0 * T**1.517083D0)
     .                * 0.980665D0 ! kgf/cm2 -> bar
        else
          He3_Pmelt = NaN
        endif
        return
      end


! T_c [mK] vs P [bar]
! Arg: P = 0 .. Pa [bar]
! Ref: Greywall. PRB33 (1986) f.5
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
        else
          He3_Tc = NaN
        endif
        return
      end

! T_ab [mK] vs P [bar]
! Arg: P = 0 .. Pmelt [bar]
! Ref: Greywall. PRB33 (1986) f.15
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
        endif
        if (P.lt.0D0.or.P.gt.He3_Pb) then
          He3_Tab = NaN
        endif
        return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Greywall -> PLTS2000 conversion T [mK] vs T [mK]
      function He3_gr2plts(t)
        implicit none
        real*8 t, he3_gr2plts
        he3_gr2plts = 0.96756D0*t + 0.031803D0
      end
      function He3_plts2gr(t)
        implicit none
        real*8 t, he3_plts2gr
        he3_plts2gr = (t - 0.031803D0)/0.96756D0
      end

! PLTS2000 melting pressure temperature scale [bars] vs T [mK]
! Arg: T = 0.0009 .. 0.250 [K]
      function He3_Pmelt_plts(T)
        implicit none
        include 'he3.fh'
        real*8 T
        if (T.gt.9D-4.and.T.le.0.25D0) then
          He3_Pmelt_plts =
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
         He3_Pmelt_plts = He3_Pmelt_plts * 10D0 ! MPa -> bar
        else if (T.gt.0.25D0.and.T.le.0.5D0) then
!         Interpolation (see tests/pmelt_interp.m)
          He3_Pmelt_plts =
     .      -117.241694D0 * T**5
     .      +319.841794D0 * T**4
     .      -304.391790D0 * T**3
     .      +167.316047D0 * T**2
     .       -49.031497D0 * T
     .       +34.882895D0
        else
          He3_Pmelt_plts = He3_Pmelt(T)
        endif
        return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! B_ab [mk] vs P, ttc
!! Inseob Hahn PhD thesis, p79
!      function He3_Bab(P,ttc)
!        implicit none
!        include 'he3.fh'
!        real*8 P,ttc
!        real*8 p0,p1,p2,q1,pp
!        real*8 f1,f2,f3,f4,f5
!        real*8 Bc,B0,gr, Pa,Ta, BB,X2
!        Pa = 34.338
!        pp = P/Pa
!        Bc=(3391D0 + 21500D0*pp - 8490D0*pp**2) / (1+2.098D0*pp)
!        f3 = 1.41D0
!        f4 = -0.29D0 -0.41D0*pp
!        B0 = 1.97 * ttc * (1D0+he3_f0a(P))
!        gr = 0.616D0 - 1.174D0* pp + 0.301 * pp**2
!        Ta = 1D0-(B/B0)**2/gr
!        BB = (B0/Bc)**2 / 4D0 * gr
!        f5 = (1D0 - f3*Ta**6 - f4*Ta**8 - (1D0-f3-f4)*Ta**2
!     .        + (1D0+2D0*f3+3D0*f4)*(Ta**4-Ta**2))
!     .       /(BB * (Ta**4-Ta**2)) - 1D0
!        f2 = BB * (1D0+f5) - (1D0+2D0*f3+3D0*f4)
!        f1 = 1D0 + f5-f2-f3-f4
!        X2 = dsqrt(
!     .   (f1*ttc**2 + f2*ttc**4 + f3*ttc**6 + f4*ttc**8)
!     .    / (1+f5*ttc**2))
!        He3_Bab = dsqrt(1D0 - X2**2) * Bc
!      end

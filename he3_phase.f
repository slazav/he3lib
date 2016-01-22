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
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Melting pressure [bars] vs T [K], original datasets
!     Greywall. PRB33 7520 (1986) f.A1,  T = 0.0009 - 0.25 K
!      see also tests/phase/pmelt_interp.m
      function He3_pmelt_greywall_org(T)
        implicit none
        real*8 T, He3_pmelt_greywall_org
          He3_pmelt_greywall_org = 34.3380D0
     .     - 0.19652970D-1 * (T*1D3)**(-3)
     .     + 0.61880268D-1 * (T*1D3)**(-2)
     .     - 0.78803055D-1 * (T*1D3)**(-1)
     .     + 0.13050600D0
     .     - 0.43519381D-1 * (T*1D3)
     .     + 0.13752791D-3 * (T*1D3)**2
     .     - 0.17180436D-6 * (T*1D3)**3
     .     - 0.22093906D-9 * (T*1D3)**4
     .     + 0.85450245D-12* (T*1D3)**5
      end
!     PLTS-2000, T = 0.0009 - 1 K
      function He3_pmelt_plts_org(T)
        implicit none
        real*8 T, He3_pmelt_plts_org
          He3_pmelt_plts_org = (
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
     .     - 4.5875709D0 * T**9 )
     .           * 10D0 ! MPa -> bar
      end
!     Osborne, Abraham, Weinstock, 1951, 1952,  0.5-1.5 K
      function He3_pmelt_osborne_org(T)
        implicit none
        real*8 T, He3_pmelt_osborne_org
        He3_pmelt_osborne_org =
     .    (26.8D0 + 13.1D0 * T**2)
     .              * 1.01325D0 ! atm -> bar
      end
!     Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955), 2-31K
      function He3_pmelt_mills_org(T)
        implicit none
        real*8 T, He3_pmelt_mills_org
        He3_pmelt_mills_org =
     .    (25.16D0 + 20.08201D0 * T**1.517083D0)
     .              * 0.980665D0 ! kgf/cm2 -> bar
      end
!     Interpolation function between PLTS, osborne, mills
!     see tests/phase/pmelt_interp.m
      function He3_pmelt_interp(T)
        implicit none
        real*8 T, He3_pmelt_interp
        He3_pmelt_interp =
     .   -2.8399D0*T**3 + 21.1585D0*T**2 - 3.5740D0*T + 25.1894D0
      end

! Smooth function for T = 0.0009 .. 31 K, PLTS-2000
      function He3_pmelt_plts(T) ! PLTS version
        implicit none
        include 'he3.fh'
        real*8 T
        real*8 He3_pmelt_plts_org
        real*8 He3_pmelt_mills_org
        real*8 He3_pmelt_interp
        if (T.ge.9D-4.and.T.le.1D0) then
          He3_pmelt_plts = He3_pmelt_plts_org(T)
        else if (T.gt.1D0.and.T.lt.1.1D0) then
          He3_pmelt_plts =
     .      ( (T-1.0D0)*He3_pmelt_interp(T)
     .      + (1.1D0-T)*He3_pmelt_plts_org(T) ) / (1.1D0-1.0D0)
        else if (T.ge.1.1D0.and.T.le.2.0D0) then
          He3_pmelt_plts = He3_pmelt_interp(T)
        else if (T.gt.2.0D0.and.T.lt.3.0D0) then
          He3_pmelt_plts =
     .      ( (T-2.0D0)*He3_pmelt_mills_org(T)
     .      + (3.0D0-T)*He3_pmelt_interp(T) ) / (2.0D0-3.0D0)
        else if (T.ge.3.0D0.and.T.le.31D0) then
          He3_pmelt_plts = He3_pmelt_mills_org(T)
        else
          He3_pmelt_plts = NaN
        endif
      end
! Smooth function for T = 0.0009 .. 31 K, Greywall-86
      function He3_pmelt(T)
        implicit none
        include 'he3.fh'
        real*8 T
        real*8 He3_pmelt_greywall_org
        if (T.ge.9D-4.and.T.le.0.25D0) then
          He3_pmelt = He3_pmelt_greywall_org(T)
        else if (T.gt.0.25D0.and.T.lt.0.27D0) then
          He3_pmelt =
     .      ( (T-0.25D0)*He3_pmelt_greywall_org(T)
     .      + (0.27D0-T)*He3_pmelt_plts(T) ) / (0.27D0-0.25D0)
        else if (T.ge.0.27D0.and.T.le.31D0) then
          He3_pmelt = He3_pmelt_plts(T)
        else
          He3_pmelt = NaN
        endif
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
! see tests/phase/gr2plts.m
!  t2 = (1D0 - 5.14809D-2 + 3.05134D-2*t1 - 6.9936D-3*t1**2)*t1
!  p2 = (1D0 + 8.2454D-05)*p1
!  c = (-1.2297e-02, -4.2958e-04, 8.2456e-05)
!
! T_c and Tab scaled from Greywall to PLTS
      function He3_tc_plts(P)
        implicit none
        include 'he3.fh'
        real*8 P,P1,T1
        P1 = P/(1D0 + 8.2454D-05)
        T1 = he3_tc(P1)
        He3_tc_plts =
     .    (1D0 - 5.14809D-2 + 3.05134D-2*T1 - 6.9936D-3*T1**2)*T1
      end
      function He3_tab_plts(P)
        implicit none
        include 'he3.fh'
        real*8 P,P1,T1
        P1 = P/(1D0 + 8.2454D-05)
        T1 = he3_tab(P1)
        He3_tab_plts =
     .    (1D0 - 5.14809D-2 + 3.05134D-2*T1 - 6.9936D-3*T1**2)*T1
      end

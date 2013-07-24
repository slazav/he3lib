! Melting pressure [bars] vs T [mK]
! Arg: T = 0 .. 31000 [mK]
! Note: Wrong for T<0.2mK?? - check range!
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520
! Ref: Osborne, Abraham, Weinstock, 1951, 1952
! Ref: Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955))
! Origin: Mukharskii, Dmitriev, Zavjalov

! see also: Johnson, Symko, Weatley -- Phis. Rev. Lett. 23, 1017 (1969)
      function He3_Pmelt(T)
        implicit none
        include 'he3.fh'
        real*8 T
        if (T.gt.0D0.and.T.le.250D0) then
!         Greywall. Phys. Rev.B v.33 #11 p.7520
          He3_Pmelt = He3_Pa
     .     - 0.19652970D-1*T**(-3)
     .     - 0.61880268D-1*T**(-2)
     .     - 0.78803055D-1*T**(-1)
     .     + 0.13050600D0
     .     - 0.43519381D-1*T
     .     + 0.13752791D-3*T**2
     .     - 0.17180436D-6*T**3
     .     - 0.22093906D-9*T**4
     .     + 0.85450245D-12*T**5
        else if (T.gt.250D0.and.T.le.500D0) then
!         Interpolation
          He3_Pmelt = 33.2505D0
     .     - 24.6026D0 *(T/1D3)
     .     + 34.5713D0 *(T/1D3)**2
     .     + 20.4086D0 *(T/1D3)**3
     .     - 26.3258D0 *(T/1D3)**4
        else if (T.gt.500D0.and.T.le.1500D0) then
!         Osborne, Abraham, Weinstock, 1951
!         Pm = 26.8 + 13.1 T^2 [atm], T = 0.5 .. 1.5
!         atm->bar: 1.01325
          He3_Pmelt = (26.8D0 + 13.1D0 * (T/1D3)**2) * 1.01325D0
        else if (T.gt.1500D0.and.T.le.2060D0) then
!         Interpolation
          He3_Pmelt = 24.787D0
     .    + 9.4456D0 *(T/1D3)
     .    + 0.0933D0 *(T/1D3)**2
     .    + 7.6155D0 *(T/1D3)**3
     .    - 1.5327D0 *(T/1D3)**4
        else if (T.gt.2060D0.and.T.le.31000D0) then
!         Mills, Grilly, 1955 (Phys. Rev. 99, 480486 (1955))
!         Pm = 25.16 + 20.08201 T^1.517083  [kg/cm2] P=76-3500
!         He4: -17.80 + 17.31457 T^1.555414 [kg/cm2] P=37-3500
!         kgf/cm2 -> bar: 0.980665
          He3_Pmelt = (25.16D0 + 20.08201D0 * (T/1D3)**1.517083D0)
     .                * 0.980665D0
        else
          He3_Pmelt = NaN
        endif
        return
      end



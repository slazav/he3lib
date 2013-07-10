! Melting pressure [bars] vs T [mK]
! Arg: T = 0 .. 250 [mK]
! Note: Wrong for T<0.2mK?? - check range!
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520
! Origin: Mukharskii, Dmitriev

      function He3_Pmelt(T)
        implicit none
        include 'he3.fh'
        real*8 T
        if (T.gt.0D0.and.T.le.250D0) then
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
        else
          He3_Pmelt = -1D0
        endif
        return
      end

! P=a+b*T**c
! a=25.16 b=20.08201 c=1.517083
! range P = 76..3500 kgs/cm2
! R. L. Mills and E. R. Grilly 
! Melting Curves of He3, He4, H2, D2, Ne, N2, and O2 up to 3500 kg/cm2
! Phys. Rev. 99, 480486 (1955)

! Osborne, Abraham, Weinstock, 1951, 1952 T = 0.2...1.5 K

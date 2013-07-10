! T_ab [mK] vs P [bar]
! Arg: P = 0 .. Pmelt [bar]
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520
! Origin: Mukharskii, Dmitriev

      function He3_Tab(P)
        implicit none
        include 'he3.fh'
        real*8 P, Pr
        Pr = P - He3_Pabn
        if (Pr.lt.0D0) then
          He3_Tab=He3_Tc(P)
        else
          He3_Tab= He3_Tabn
     .       - .10322623D-1*Pr
     .       - .53633181D-2*Pr**2
     .       + .83437032D-3*Pr**3
     .       - .61709783D-4*Pr**4
     .       + .17038992D-5*Pr**5
        endif
        if (P.lt.0D0.or.P.gt.He3_Pmelt(He3_Tab)) then
          He3_Tab = -1D0
        endif
        return
      end

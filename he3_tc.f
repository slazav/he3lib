! T_c [mK] vs P [bar]
! Arg: P = 0 .. Pa [bar]
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520
! Origin: Mukharskii, Dmitriev

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

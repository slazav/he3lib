! Molar Volume [cm**3/mole] vs P [bar]
! Arg: P = 0 .. Pa [bars]
! Note: No temperature dependance,
!       no check for solid phase.
! Ref: Greywall. Ph.Rev.B v.33 #11 p.7520 ref. 27
!      That is Wheatley Rev.Mod.Phys. 47,415(1975)
! Origin: Mukharskii, Dmitriev

      function He3_Vm(P)
        implicit none
        include 'he3.fh'
        real*8 P
        if (P.gt.0D0.and.P.le.He3_Pa) then
          He3_Vm=   36.837231D0
     .       - 1.1803474D0*P
     .       + 0.0834214D0*P**2
     .       - 0.3885962D-2*P**3
     .       + 0.94759780D-4*P**4
     .       - 0.91253577D-6*P**5
        else
          He3_Vm = -1D0
        endif
      end

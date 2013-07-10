! R-Gas constant GAMMA=C/RT [1/(K*mol)] vs P [bar]
! Greywall. Phys. Rev.B v.33 #11 p.7520
! Origin: Mukharskii, Dmitriev

      function He3_gammaf(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_gammaf = .27840464D+1
     .            + .69575243D-1*P
     .            - .14738303D-2*P**2
     .            + .46153498D-4*P**3
     .            - .53785385D-6*P**4
      end

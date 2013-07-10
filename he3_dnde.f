! Attention density of state must be multiplyied by ANA later.
! Origin: Mukharskii, Dmitriev

      function He3_dNdE(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_dNdE=3D0*He3_gammaf(P)/AKB/PI**2
      end

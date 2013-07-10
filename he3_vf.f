! Fermi velocity [cm/s] vs P [bar]
! Origin: Mukharskii, Dmitriev

      function He3_Vf(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_Vf=He3_Pf(P)/He3_Meff(P)
      end

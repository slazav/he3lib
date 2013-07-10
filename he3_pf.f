! Fermi momentum [sgs] vs P [bar]
! Origin: Mukharskii, Dmitriev

      function He3_Pf(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_Pf=HC*(3D0*PI**2*ANA/He3_Vm(P))**.3333333D0
      end

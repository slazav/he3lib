! Attention density of state must be multiplyied by ANA later.
! Origin: Mukharskii, Dmitriev

      function He3_dNdE(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_dNdE=3D0*He3_gammaf(P)/AKB/PI**2
      end

! Effective mass [g] vs P [bar].
! Origin: Mukharskii, Dmitriev

      function He3_Meff(P)
        implicit none
        include 'he3.fh'
        real*8 P,PF
        PF=He3_Pf(P)
C       He3_Meff=He3_gammaf(P)*R*HC*He3_Vm(P)*(HC/PF)*3
C       print *,He3_gammaf(P),He3_gammaf(P)*R*HC,R
        He3_Meff=He3_dNdE(P)*PF/3D0*PF
C       He3_Meff=He3_dNdE(P)*PF/3./He3_Vm(P)*PF
      end

! Fermi momentum [sgs] vs P [bar]
! Origin: Mukharskii, Dmitriev

      function He3_Pf(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_Pf=HC*(3D0*PI**2*ANA/He3_Vm(P))**.3333333D0
      end

! Fermi velocity [cm/s] vs P [bar]
! Origin: Mukharskii, Dmitriev

      function He3_Vf(P)
        implicit none
        include 'he3.fh'
        real*8 P
        He3_Vf=He3_Pf(P)/He3_Meff(P)
      end

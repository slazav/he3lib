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

! Suseptibility [sgs] vs P [bar], T [mK]
! Origin: Mukharskii, Dmitriev

      function He3_susept(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T,Y,TTC,Z0
        He3_susept = he3_chi0(P)
        TTC=T/He3_Tc(P)
        if (TTC.LT.1D0) then
          Z0 = He3_z0(P)
          Y  = He3_yosida(TTC)
C         print *,Z,Y,((1.+Z0/4.)*(2.+Y)/3.)/(1.+Z0/4.*(2.+Y)/3.)
          He3_susept = He3_susept *
     .      ((1D0+Z0/4D0)*(2D0+Y)/3D0)/(1D0+Z0/4D0*(2D0+Y)/3D0)
        end if
      end

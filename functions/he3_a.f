!HH> A phase

!> Legget frequency nu_a [Hz] vs P, ttc
!> Interpolation formula by A.Yudin based on Halperin and ROTA data
!> There is no check that t > t_ab
      function he3_nu_a(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p

        if (ttc.gt.1D0) then
          he3_nu_a = 0D0
        else
          he3_nu_a = (0.37155D0 + 0.15366D0*p -0.00113D0* P**2)
     .     * ((1D0-ttc) - 1.14915D0*(1-ttc)**2 +0.43495D0*(1-ttc)**3)
          he3_nu_a = sqrt(he3_nu_a*1D10)
        endif
      end

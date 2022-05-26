!HH> 3He-4He mixtures

!H>  Constants
      block data he34_const_block
        implicit none
        include 'he3.fh'

        data
     .    he34_xcr  /0.674D0/, !C> critical point at saturated pressure, concentration of He3
     .    he34_tcr  /0.867D0/  !C> critical point at saturated pressure, temperature [K]
      end


!> Phase diagram functions are from Chaudhry PhD thesis (Massachusetts, 2009)

!H> Phase diagram at saturated vapor pressure

!> Dilute phase-separation curve, concentration vs temperature [K]
!> (0.08014 < x < 0.674; 0.15 K < T < 0.867 K):
!> Chaudhry PhD thesis (Massachusetts, 2009)
!> Below 0.15 K: Edwards, Ifft, Sarwinski, Phys.Rev. 177, 380 (1969) Eq.23
      function he34_xdil(t) !F>
        implicit none
        include 'he3.fh'
        real*8 t, dt
        if (t.gt.he34_tcr) then
          he34_xdil = NaN
        elseif (t.gt.0.15) then
          dt = t-he34_tcr
          he34_xdil = he34_xcr
     .      - 0.209148D0*dt/(dt-0.080280D0)
     .      + 0.960222D0*dt + 0.549920D0*dt**2
        else
!          ! formula from Lancaster wire programs
!          he34_xdil = 0.066D0 + 0.5056D0*t**2
!     .      - 0.2488D0*t**3 + 18.22D0*t**4 - 74.22D0*t**5
          ! Edwards-1969, eq.23
          he34_xdil = 0.064D0*(1D0 + 10.8D0*t**2)
        endif
      end

!> Concentrated phase-separation curve, concentration vs temperature [K]
!> (0.674 < x < 1; 0.15 K < T < 0.867 K).
!> Chaudhry PhD thesis (Massachusetts, 2009)
!> Below 0.15 K: Edwards, Ifft, Sarwinski, Phys.Rev. 177, 380 (1969) Eq.25.
      function he34_xcon(t) !F>
        implicit none
        include 'he3.fh'
        real*8 t, dt
        if (t.gt.he34_tcr) then
          he34_xcon = NaN
        elseif (t.gt.0.15) then
          dt = t-he34_tcr
          he34_xcon = he34_xcr + 0.316170D0*dt**3
     .       - 0.180743D0*dt**2 - 0.746805D0*dt
        else
          ! Edwards-1969, eq.25
          he34_xcon = 1D0 - 1.13D0*t**1.5D0*dexp(-0.71D0/t)
        endif
      end

!> Lambda curve, temperature [K] vs concentration
!> Chaudhry PhD thesis (Massachusetts, 2009)
      function he34_tlambda(x) !F>
        implicit none
        include 'he3.fh'
        real*8 x, dx
        if (x.gt.he34_xcr.or.x.lt.0D0) then
          he34_tlambda = NaN
        else
          dx = x - he34_xcr
          he34_tlambda = he34_tcr
     .      - 2.620259D0*dx - 1.023726D0*dx**2
        endif
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!H> Phase diagram at pressures 0..10bar

!> Tricritical line vs pressure [bar], 0..10 bar
!> Chaudhry PhD thesis (Massachusetts, 2009)
      function he34_tcr_p(p) !F>
        implicit none
        include 'he3.fh'
        real*8 p
        if (p.lt.0D0.or.p.gt.10D0) then
          he34_tcr_p = NaN
        else
          he34_tcr_p = he34_tcr
     .      - 0.12992576D0*p/(p+2.5967345D0) - 6.457263D-4*p
        endif
      end

!> Concentration along the tricritical line vs pressure [bar]:
!> Chaudhry PhD thesis (Massachusetts, 2009)
      function he34_xcr_p(p) !F>
        implicit none
        include 'he3.fh'
        real*8 p, dt
        if (p.lt.0D0.or.p.gt.10D0) then
          he34_xcr_p = NaN
        else
          dt = he34_tcr - he34_tcr_p(p)
          he34_xcr_p = he34_xcr
     .     + 0.3037124D0*dt - 4.41225D6*dt**9
        endif
      end

!> Dilute phase-separation curve, x vs t[K] and p[bar]
!> Chaudhry PhD thesis (Massachusetts, 2009)
      function he34_xdil_p(t,p) !F>
        implicit none
        include 'he3.fh'
        real*8 t, p, dt, K0,K1,K2,Ka
        if (t.lt.0.15D0.or.t.gt.he34_tcr_p(p)) then
          he34_xdil_p = NaN
        elseif (p.lt.0D0.or.p.gt.10D0) then
          he34_xdil_p = NaN
        else
          dt = t-he34_tcr_p(p)
          K0 = -0.209148D0 - 0.1269791D0*p + 0.0102283D0*p**2
          K1 =  0.960222D0 - 0.2165742D0*p + 0.0169801D0*p**2
          K2 =  0.549920D0 - 0.1198491D0*p + 0.0092997D0*p**2
          Ka =  0.080280D0 +0.02291499D0*p - 0.0020886D0*p**2
          he34_xdil_p = he34_xcr_p(p)
     .     + K0*dt/(dt-Ka) + K1*dt + K2*dt**2
        endif
      end

!> Concentrated phase-separation curve, x vs t[K] and p[bar]
!> Chaudhry PhD thesis (Massachusetts, 2009)
      function he34_xcon_p(t,p) !F>
        implicit none
        include 'he3.fh'
        real*8 t, p, dt, K1,K2,K3
        if (t.lt.0.15D0.or.t.gt.he34_tcr_p(p)) then
          he34_xcon_p = NaN
        elseif (p.lt.0D0.or.p.gt.10D0) then
          he34_xcon_p = NaN
        else
          dt = t-he34_tcr_p(p)
          K1 = -0.746805D0 + 0.0173549D0*p - 0.0028598D0*p**2
          K2 = -0.180743D0 + 0.1120251D0*p - 0.0152076D0*p**2
          K3 =  0.316170D0 + 0.1723264D0*p - 0.0201411D0*p**2
          he34_xcon_p = he34_xcr_p(p)
     .     + K1*dt + K2*dt**2 + K3*dt**3
        endif
      end

!> Lambda curve, t[K] vs x and p[bar]
!> Chaudhry PhD thesis (Massachusetts, 2009)
      function he34_tlambda_p(x,p) !F>
        implicit none
        include 'he3.fh'
        real*8 x, p, dx, K1,K2
        if (x.gt.he34_xcr_p(p).or.x.lt.0D0) then
          he34_tlambda_p = NaN
        elseif (p.lt.0D0.or.p.gt.10D0) then
          he34_tlambda_p = NaN
        else
          dx = x - he34_xcr_p(p)
          K1 = -2.620259D0 + 0.0024823D0*p + 0.0009255D0*p**2
          K2 = -1.023726D0 + 0.0013175D0*p + 0.0009397D0*p**2
          he34_tlambda_p = he34_tcr_p(p) + K1*dx + K2*dx**2
        endif
      end


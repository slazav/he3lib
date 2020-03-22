!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Helium-4 parameters

      block data he4_const_block
        implicit none
        include 'he3.fh'

        data
     .    he4_amass  /6.6465D-24/,   ! He4 atom mass, [g]
     .    he4_mmass  /4.002602D0/    ! He4 molar mass, [g/mol]
!         Swenson-1952 temperature data is shifted to have Tc(Pvap) = 2.1720K (1958 temperature scale)
!         in 1958 temperature scale umHg = 1.33322387e-6 bar
     .    he4_tcv    /2.1720D0/,     ! superfluid transition temperature at vapor pressure [K] (1958 temperature scale)
     .    he4_pcv    /0.050396D0/,   ! vapor pressure at superfluid transition [bar] (1958 temperature scale)
     .    he4_tcm    /1.750D0/,      ! superfluid transition temperature at melting pressure [K] (Swenson-1952, shifted)
     .    he4_pcm    /30.033D0/,     ! superfluid transition pressure at melting curve [bar] (Swenson-1952)
     .    he4_tcr    /5.1994D0/,     ! critical temperature [K] (1958 temperature scale)
     .    he4_pcr    /2.2905D0/      ! critical pressure [bar] (1958 temperature scale)
      end

! Superfluid transition temperature [K] vs pressure [bar]
! Swenson-1952, table1 (1 atm = 1.01325 bar)
! Values are shifted to get Tc = 2.1720 at vapor pressure (as in 1958 temperature scale)
      function he4_tc(p)
        implicit none
        include 'he3.fh'
        real*8 p

        if (p.lt.he4_pcv.or.p.gt.he4_pcm) then
          he4_tc = NaN
        else
          he4_tc = 2.1862228445D0
     .           - 0.0109123116D0 * p
     .           - 0.0001037245 * p**2
     .           - 0.014  ! shift: 2.1720 - 2.186
        endif
      end

! Melting pressure [bar] vs temperature [K], 0..4K
! Swenson-1950,1951
      function he4_pmelt(t)
        implicit none
        include 'he3.fh'
        real*8 t

        if (t.eq.0D0) then
          he4_pmelt = 25.00D0
          return
        endif
        if (t.gt.0D0.and.t.le.4D0) then
          he4_pmelt = 25.00D0 + 0.33046D0
     .       * dexp(15.901D0/t - 24.1085D0/t**2 + 0.820147D0*t)
          return
        endif
        he4_pmelt = NaN
      end

! Vapor pressure [bar] vs temperature [K], 0..tcr
! Fit of 1958 temperature scale (~1.5% accuracy)
      function he4_pvap(t)
        implicit none
        include 'he3.fh'
        real*8 t
        real*8 a,b,c,d,e,ee,k1,k2

        if (t.eq.0D0) then
          he4_pvap = 0D0
          return
        endif

        if (t.le.he4_tcr) then

          a = 40.896039500269D0
          b = 6.73825999245785D0
          c = -3.147638667442D0
          d = 3.97695025409929D0
          e = -1.77535194090704D0
          ee = 3.08578394156274D0
          k1 = 0.0606579218927369D0
          k2 = 3.42982579814079D0
          he4_pvap = a*(t/b)**k2 * dexp(-b/t)
     .       * (1D0 + c*(t/b) + d*(t/b)**2 + e*(t/b)**3)
     .       * (1D0 + ee*(t/b)**k1)

!          ! old fit version, 3% acc
!          he4_pvap = dexp(-2.38491D0 - 8.14649D0/t
!     .       + 2.19785D0*t - 0.423192D0*t**2 + 0.0342826D0*t**3)
          return
        endif

        he4_pvap = NaN
      end

! Molar volume [cm^3/mol] vs T [K], saturated vapor pressure
! Kerr and Taylor 1964
      function he4_vm(t)
        implicit none
        include 'he3.fh'
        real*8 t, dt

!       2.22K .. 4.4K
!       In the paper there is another expression which uses he3_pvap
!       This fit is also good but has some small derivative step
        dt = t - 2.22D0
        if (dt.gt.0D0) then
          he4_vm = 10D0**(1.43771882061654D0
     .           + 0.00861848D0*dt
     .           + 0.0127383D0*dt**2)
          return
        endif

!       Tc .. 2.22K
!       fitting formula from the paper
        dt = t - 2.1720D0 ! he4_tcv
        if (dt.gt.0D0.and.t.le.2.22D0) then
          he4_vm = 10D0**(1.4375451D0
     .           + 0.013287D0*dt
     .           + 0.007331D0*dt*dlog(dt)/dlog(10D0))
          return
        endif

!       Tc
        if (dt.eq.0D0) then
          he4_vm = 10D0**(1.4375451D0)
          return
        endif

!       2.15K .. Tc
!       fitting formula from the paper
        if (dt.lt.0D0.and.t.ge.2.15D0) then
          he4_vm = 10D0**(1.4375451D0
     .           - 0.002102D0*dt
     .           + 0.007313D0*dt*dlog(-dt)/dlog(10D0))
          return
        endif

!       0 .. 2.15 K
!       This fit is good, but has some small step at 2.15
        if (t.gt.0.and.t.lt.2.15D0) then
          he4_vm = 27.5793D0
     .           - 55.912D0 * exp(-9.8751D0/t)
     .           + 0.598036D0 * exp(-(t-2.67447D0)**2D0/0.659719D0)
          return
        endif

!       0
        if (t.gt.0.and.t.lt.2.15D0) then
          he4_vm = 27.5793D0
          return
        endif 
        he4_vm = NaN
      end



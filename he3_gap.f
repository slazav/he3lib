! BCS gap / (kB Tc) for pure 3He-B, t = T / Tc
! Newton iteration based on a note by EVT & RH
! From ROTA texture library
      function he3_bcsgap(ttc)
        implicit none
        include 'he3.fh'
        integer n, m
        real*8 ttc,root,y,dy,ynew,g,dg
        if (ttc.ge.1) then
          he3_bcsgap = 0D0
          return
        end if
        if (ttc.lt.0) then
          he3_bcsgap = NaN
          return
        end if
        m = 30
        dy = 1.0D0
        ynew = 1.7638D0*SQRT(1D0-ttc)/(2D0*const_pi)
        do while (ABS(dy) > 1.0D-8)
          y = ynew
          root = SQRT((dble(m)*ttc)**2+y**2)
          g = LOG((dble(m)*ttc+root) / (2D0*dble(m)))
     .        - (1D0/dble(m)**2 - dble(m)*(ttc/root)**3)/24D0
          dg = y/(root*(dble(m)*ttc+root))
     .        - dble(m)*ttc**3*y/(8D0*root**5)
          do n=1,m
            root=SQRT((ttc*(dble(n)-0.5D0))**2 + y**2)
            g = g + 1D0/(dble(n)-0.5D0) - ttc/root
            dg = dg + ttc*y/root**3
          end do
          dy = g/dg
          ynew = ynew-dy
        end do
        he3_bcsgap = 2D0*const_pi*ynew
        if (ttc.ge.1D0) he3_bcsgap = 0D0
      end

! BCS gap / (kB Tc) for pure 3He-B, t = T / Tc
! Einzel approximation (D.Einzel JLTP 84 (1991))
! <0.5% accuracy in the whole temperature range
! 70 times faster
      function he3_bcsgap_fast(ttc)
        implicit none
        include 'he3.fh'
        real*8 ttc, dsc, ccn, c1,c2
        if (ttc.ge.1) then
          he3_bcsgap_fast = 0D0
          return
        end if
        if (ttc.lt.0) then
          he3_bcsgap_fast = NaN
          return
        end if
        he3_bcsgap_fast = 1.764D0 *
     .    dtanh( const_pi/1.764D0 *
     .      sqrt(2D0/3D0 * 1.426D0 * (1D0-ttc)/ttc *
     .          (1+0.1916D0*(1-ttc) + 0.2065D0*(1-ttc)**2)))
      end

! Trivial strong-coupling correction to the BCS energy gap
      function he3_trivgap(ttc,p)
        implicit none
        include 'he3.fh'
        integer it, ic
        real*8 ttc, p, dcpcn, wt1, wt2, wc1, wc2, corr
        real*8 c,x
        dimension c(5), x(11,5)
        c = (/ 1.43D0,1.6D0,1.8D0,2.D0,2.2D0 /)
        dcpcn = 41.9D0 / he3_vm(p) + 0.322D0
        x(11,1:5) = (/ 1.D0,1.056D0,1.115D0,1.171D0,1.221D0 /)
        x(10,1:5) = (/ 1D0,1.048D0,1.097D0,1.141D0,1.18D0 /)
        x(9,1:5)  = (/ 1D0,1.041D0,1.083D0,1.119D0,1.15D0 /)
        x(8,1:5)  = (/ 1D0,1.036D0,1.072D0,1.102D0,1.128D0 /)
        x(7,1:5)  = (/ 1D0,1.032D0,1.063D0,1.089D0,1.112D0 /)
        x(6,1:5)  = (/ 1D0,1.028D0,1.056D0,1.079D0,1.099D0 /)
        x(5,1:5)  = (/ 1D0,1.026D0,1.051D0,1.073D0,1.091D0 /)
        x(4,1:5)  = (/ 1D0,1.024D0,1.049D0,1.069D0,1.086D0 /)
        x(3,1:5)  = (/ 1D0,1.024D0,1.048D0,1.068D0,1.085D0 /)
        x(2,1:5)  = (/ 1D0,1.024D0,1.048D0,1.068D0,1.085D0 /)
        x(1,1:5)  = (/ 1D0,1.024D0,1.048D0,1.068D0,1.085D0 /)
        it=INT(ttc*10D0 - 1D-5) + 1
        wt1 = (0.1D0*dble(it)-ttc)/0.1D0
        wt2 = 1D0 - wt1
        ic = 1
        do
          if (dcpcn < c(ic+1) ) exit
          ic = ic+1
          if (ic == 4) exit
        end do
        wc1 = (c(ic+1)-dcpcn)/(c(ic+1)-c(ic))
        wc2 = 1D0 - wc1
        corr = wt1*(wc1*x(it,ic)+wc2*x(it,ic+1))
        corr = corr + wt2*(wc1*x(it+1,ic)+wc2*x(it+1,ic+1))
        he3_trivgap = he3_bcsgap(ttc)*corr
      end

      function he3_z3(ttc,gap)
        implicit none
        include 'he3.fh'
        integer i,maxi
        real*8 ttc,gap,y,help,mt,corr1,corr2,sum
        y = gap/(2D0*const_pi)
        sum = 0D0
        maxi = 100
        mt = dble(maxi)*ttc
        do i=1,maxi
           sum = sum + ttc/((ttc*(dble(i)-0.5D0))**2+y**2)**1.5D0
        end do
        help = SQRT(mt**2 + y**2)
        corr1 = 1D0/(help*(mt+help))
        corr2 = mt**3/(8D0*help**5)
        he3_z3 = y**2*(sum + corr1 - corr2)
      end

      function he3_z5(ttc,gap)
        implicit none
        include 'he3.fh'
        integer i,maxi
        real*8 ttc,gap,y,help,mt,corr1,corr2,sum
        y = gap/(2D0*const_pi)
        sum = 0D0
        maxi = 100
        mt = dble(maxi)*ttc
        do i=1,maxi
           sum = sum + ttc/((ttc*(dble(i)-0.5D0))**2 + y**2)**2.5D0
        end do
        help = SQRT(mt**2+y**2)
        corr1 = (mt+2D0*help)/(3D0*help**3*(mt+help)**2)
        corr2 = 5D0*mt**3/(24D0*help**7)
        he3_z5 = y**4*(sum + corr1 - corr2)
      end

      function he3_z7(ttc,gap)
        implicit none
        include 'he3.fh'
        integer i,maxi
        real*8 ttc,gap,y,help,mt,corr1,corr2,sum
        y = gap/(2D0*const_pi)
        sum = 0D0
        maxi = 100
        mt=dble(maxi)*ttc
        do i=1,maxi
           sum = sum + ttc/((ttc*(dble(i)-0.5D0))**2+y**2)**3.5D0
        end do
        help = SQRT(mt**2+y**2)
        corr1 = (11D0*mt*mt + 9D0*mt*help + 8D0*y*y)/
     .          (15D0*help**5*(mt+help)**3)
        corr2 = 7D0*mt**3/(24D0*help**9)
        he3_z7 = y**6*(sum + corr1 - corr2)
      end


! Yosida function vs T/Tc, gap -- old code
      function he3_yosida0(ttc,gap)
        implicit none
        include 'he3.fh'
        real*8 ttc, gap
        he3_yosida0 = 1D0 - he3_z3(ttc,gap)
      end


! Integrand for Yosida function calculations
! x = tanh(\xi)/2T change is made to get good integrand
! and [0:1] integrating range.
      function he3_yosida_int(x, ttc,gap, n)
        implicit none
        include 'he3.fh'
        real*8 x, ttc,gap, xi, ek, n
        real*8 he3_yosida_int
        xi = datanh(x)*2D0*ttc
        ek=dsqrt(xi**2 + gap**2)
        he3_yosida_int =
     .     (xi/ek)**n
     .   / (dcosh(ek/(2D0*ttc)))**2
     .   * dcosh(xi/(2D0*ttc))**2
     .   * 2D0*ttc
      end

! Yosida function vs T/Tc, gap
! See D.Einzel JLTP 84
      function he3_yosida(ttc, gap, n)
        implicit none
        include 'he3.fh'
        real*8 ttc, gap, n
        real*8 dx, xp, xm
        real*8 he3_yosida_int
        integer i, maxi
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_yosida=NaN
          return
        endif
        maxi=100
        dx=1D0/dble(maxi)
        he3_yosida = 0D0
        ! intergation of he3_yosida_int from 0 to 1 using Gaussian quadrature
        do i=1,maxi 
          xp = dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          he3_yosida = he3_yosida
     .       + he3_yosida_int(xp, ttc, gap, n) * dx/2D0
     .       + he3_yosida_int(xm, ttc, gap, n) * dx/2D0
        enddo
        he3_yosida = he3_yosida / (2D0*ttc)
      end

!! Yosida0 -- does not work
!! Einzel approximation (D.Einzel JLTP 130 (2003))
!      function he3_yosida0_fast(ttc, gap)
!        implicit none
!        include 'he3.fh'
!        real*8 ttc, gap, gap0, k
!        if (ttc.lt.0D0.or.ttc.gt.1D0) then
!          he3_yosida0_fast=NaN
!          return
!        endif
!        ! gap at ttc=0. scale from bcsgap to our gap function
!        gap0 = gap * he3_bcsgap_fast(0D0)/he3_bcsgap_fast(ttc)
!        k = (2.5D0 - gap0)
!     .    / (1D0 - dsqrt(2D0*const_pi*gap)
!     .             *dexp(-gap)*(1D0 + 3D0/8D0/gap))
!
!        he3_yosida0_fast =
!     .     dsqrt(2D0*const_pi*gap/ttc) * dexp(-gap/ttc)
!     .     * (1D0 + 3D0/8D0 * ttc/gap)
!     .     * (1D0 - ttc**k)
!     .   + dexp(gap - gap/ttc) * ttc**(k-0.5D0)
!      end

! Collision integral in Einzel approximation
! Einzel, Wolfle, Hirschfeld, JLTP80 (1990), Appendix, p.66
      function he3_coll_int(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 he3_coll_int, xi, ttc, gap, x
        real*8 J0,J1,J2,J3, K0,K1,K2,K3, I0,I1,I2,I3
        real*8 a0,a1,a2,a3, b0,b1,b2, g0,d0

        x = gap/ttc
        a0 = -0.5768D0
        a1 =  0.2694D0
        a2 =  0.29D0
        a3 = -0.08D0
        J0 = 3D0/4D0/const_pi**0.5D0
     .       * (1D0+2D0*x)**1.5D0
     .       / (dexp(x) + a0 + a1*x + a2*x**2 + a3*x**3)
        J1 = x**2 / (2D0*const_pi)**0.5D0
     .       * (0.5D0 + x)**1.5D0
     .       / (1.3D0 + x**2) * dexp(-x)
        J2 = x**2 / (2D0*const_pi)**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (dexp(x)+2.3D0)
        J3 = 3D0*x**4 / (2D0*const_pi)**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (0.2D0 + x**2)
     .       / (dexp(x)+1D0+0.1D0*x**2)
        b0 =  3.4296D0
        b1 = -3.2148D0
        b2 =  2.375D0
        K0 = 9D0/8D0/(2D0*const_pi)**0.5D0
     .       / (0.5D0 + x)**0.5D0
     .       / (dexp(x) + b0 + b1*x + b2*x**2)
        K1 = - 5D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**2 / (1D0+x)**0.5D0
     .       * dexp(-x) / (const_pi**2 + x**2)
        K2 = 3D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**2 / (1D0+x)**0.5D0
     .       * dexp(-x) / (127D0/150D0 * const_pi**2 + x**2)
        K3 = - 15D0/8D0/(2D0*const_pi)**0.5D0
     .       * x**4/(1D0+x**2)**0.5D0
     .       * dexp(-x) / (const_pi**2 + x**2)**2
        I0 = J0 + K0 * (xi/ttc)**2
        I1 = J1 + K1 * (xi/ttc)**2
        I2 = J2 + K2 * (xi/ttc)**2
        I3 = J3 + K3 * (xi/ttc)**2
        he3_coll_int = ( I0 - g0*(I1+I2) + d0*I3 )
      end

! Collision integral for low temp (good for < 0.7Tc)
! Einzel, JLTP84 (1991), p.345
      function he3_coll_int_lt(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 he3_coll_int_lt, xi, ttc, gap, x, g0,d0,w0
        x=xi/dsqrt(2*ttc*gap)
        w0 = (1-2D0/3D0*g0 + d0)
        he3_coll_int_lt =
     .    3D0/2D0/const_pi * gap/ttc * he3_yosida(ttc,gap, 0D0)
     .    * (w0 + ttc/gap*(0.75D0*(1D0 + x**2) * w0
     .                      - (1D0+2D0*x**2)*(g0/3D0+d0) ))
      end

! Collision integral for high temp (not very useful)
! Einzel, JLTP84 (1991), p.345
      function he3_coll_int_ht(xi,ttc, gap, g0, d0)
        implicit none
        include 'he3.fh'
        real*8 he3_coll_int_ht, xi, ttc, gap, x, g0,d0,w0
        he3_coll_int_ht = 1D0 + (xi**2 + gap**2)/(ttc*const_pi)**2
      end

! Integrand for tau_av calculations
! e = tanh(\xi)/2ttc change is made to get good integrand (e = 0..1)
      function he3_tau_av_int(x,ttc, gap, g0, d0)
        implicit none
        real*8 x,ttc, gap, g0, d0, xi, ek
        real*8 he3_coll_int, he3_tau_av_int
        xi = datanh(x)*2D0*ttc
        ek=dsqrt(xi**2 + gap**2)
        he3_tau_av_int = he3_coll_int(xi,ttc, gap, g0, d0)
     .   / (dcosh(ek/(2D0*ttc)))**2
     .   * dcosh(xi/(2D0*ttc))**2
     .   * 2D0*ttc
      end

! Quasiparticle lifetime at fermi level
! Einzel JLTP84 (1991) p.344
      function he3_tau0(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, g0, d0, tn
        real*8 he3_coll_int
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_tau0 = NaN
          return
        endif
        g0  = he3_scatt_g0(p)
        d0  = he3_scatt_d0(p)
        gap = he3_trivgap(ttc, p)
        tn  = he3_tau_n0tc(p) / ttc**2
        he3_tau0 = tn / he3_coll_int(0D0, ttc, gap, g0, d0);
      end

! Averaged quasiparticle lifetime
! Einzel JLTP84 (1991) p.345
      function he3_tau_av(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap, sum, g0, d0, Y0, tn
        real*8 dx, xp, xm
        real*8 he3_tau_av_int
        integer i, maxi
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_tau_av=NaN
          return
        endif
        g0  = he3_scatt_g0(p)
        d0  = he3_scatt_d0(p)
        gap = he3_trivgap(ttc, p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        tn  = he3_tau_n0tc(p) / ttc**2
        sum = 0D0
        maxi=1000
        dx=1D0/dble(maxi)
        ! intergation from 0 to 1 using Gaussian quadrature
        do i=1,maxi 
          xp = dx * (dble(i) - 0.5D0 + 0.5D0/dsqrt(3D0))
          xm = dx * (dble(i) - 0.5D0 - 0.5D0/dsqrt(3D0))
          sum = sum
     .       + he3_tau_av_int(xp, ttc, gap, g0, d0) * dx/2D0
     .       + he3_tau_av_int(xm, ttc, gap, g0, d0) * dx/2D0
        enddo
        he3_tau_av = Y0 * tn / sum
      end

! Spin diffusion transport time, s
! Einzel JLTP84 (1991) p.349
      function he3_tau_dperp(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, l1a, gap, y0, y2
        l1a = he3_scatt_l1a(p)
        gap = he3_trivgap(ttc, p)
        y0  = he3_yosida(ttc, gap, 0D0)
        y2  = he3_yosida(ttc, gap, 2D0)
        he3_tau_dperp = he3_tau_av(ttc,p)
     .   /(1D0 - l1a*(4D0/5D0*y0 + 1D0/5D0*y2)/y0)
      end

! He3-B suseptibility
! see VW book ch.10 p.449; ch2 p.90
! see Wheatley-75 f 3.7
! There is also additional term to 3*chi0
!  + 2/5 F2a (1-Y0)^2
      function He3_chi_b(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,G,Y,F
        F = He3_f0a(p)
        He3_chi_b = he3_chi_n(p)*(1D0+F)
        if (ttc.LT.1D0) then
          G  = he3_bcsgap(ttc)
          Y  = He3_yosida(ttc, G, 0D0)
          He3_chi_b = He3_chi_b *
     .      ((1D0+F) * (2D0+Y)/3D0) /
     .      ( 1D0+F * (2D0+Y)/3D0)
        end if
      end

! B-phase Leggett freq
      function he3_nu_b(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap
        gap  = he3_trivgap(ttc,p) * const_kb * he3_tc(p)/1D3 ! mk->K
        he3_nu_b = dsqrt(3D0 / 8D0 / const_pi / he3_chi_b(ttc,p))
     .    * he3_gyro**2 * const_hbar * he3_2n0(p) / 4D0
     .    * gap * dlog(he3_tfeff(p)*const_kB/gap)
      end

! Suseptibility [sgs] vs P [bar], T [mK]
! Origin: Mukharskii, Dmitriev

      function He3_susept(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T,G,Y,TTC,F
        He3_susept = he3_chi_n(P)
        TTC=T/He3_Tc(P)
        if (TTC.LT.1D0) then
          F = He3_f0a(P)
          G  = he3_bcsgap(ttc)
          Y  = He3_yosida0(ttc, G)
          He3_susept = He3_susept *
     .      ((1D0+F) * (2D0+Y)/3D0) /
     .      ( 1D0+F * (2D0+Y)/3D0)
        end if
      end



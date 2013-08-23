! BCS gap / (kB Tc) for pure 3He-B, t = T / Tc
! Newton iteration based on a note by EVT & RH
! From ROTA texture library
      function he3_bcsgap(ttc)
        implicit none
        include 'he3.fh'
        integer n, m
        real*8 ttc,root,y,dy,ynew,g,dg
        if (ttc.ge.1D0) then
          he3_bcsgap = 0D0
          return
        end if
        if (ttc.lt.0D0) then
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
! Einzel approximation (D.Einzel JLTP 84 (1991) f.68)
! <0.5% accuracy in the whole temperature range
! 70 times faster
      function he3_bcsgap_fast(ttc)
        implicit none
        include 'he3.fh'
        real*8 ttc, dsc, ccn, c1,c2
        if (ttc.ge.1D0) then
          he3_bcsgap_fast = 0D0
          return
        end if
        if (ttc.lt.0D0) then
          he3_bcsgap_fast = NaN
          return
        end if
        he3_bcsgap_fast = 1.764D0 *
     .    dtanh( const_pi/1.764D0 *
     .      sqrt(2D0/3D0 * 1.426D0 * (1D0-ttc)/ttc *
     .          (1D0+0.1916D0*(1D0-ttc) + 0.2065D0*(1D0-ttc)**2)))
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
! x = tanh(\xi)/2 change is made to get good integrand
! and [0:1] integrating range.  d\xi -> 2 dx / (1-x**2)
! See also tests/plot_yosida_int.m
      function he3_yosida_int(x, args)
        implicit none
        include 'he3.fh'
        real*8 x, args(3), ttc, gap, xi, ek, n, C
        real*8 he3_yosida_int
        ttc=args(1)
        gap=args(2)
        n=args(3)
        C=2D0
        xi = datanh(x)*C
        ek=dsqrt(xi**2 + gap**2)
        he3_yosida_int =
     .     (xi/ek)**n
     .     / (dcosh(ek/(2D0*ttc)))**2 / 2D0/ttc
     .     * C / (1D0-x**2)
      end


! Yosida function vs T/Tc, gap
! See D.Einzel JLTP 84
      function he3_yosida(ttc, gap, n)
        implicit none
        include 'he3.fh'
        real*8 ttc, gap, n, args(3)
        real*8 he3_yosida_int
        external he3_yosida_int

        if (ttc.lt.0D0) then
          he3_yosida=NaN
          return
        endif
        if (ttc.eq.0D0) then
          he3_yosida=0D0
          return
        endif
        if (ttc.gt.1D0) then
          he3_yosida=1D0
          return
        endif

        args = (/ttc, gap, n/)
        he3_yosida = math_dint(he3_yosida_int, 0D0, 1D0, 100, args)
      end

! Yosida0 -- does not work
! Einzel approximation (D.Einzel JLTP 130 (2003))
      function he3_yosida0_fast(ttc, gap)
        implicit none
        include 'he3.fh'
        real*8 ttc, gap, gap0, k
        if (ttc.lt.0D0.or.ttc.gt.1D0) then
          he3_yosida0_fast=NaN
          return
        endif
        ! gap at ttc=0. scale from bcsgap to our gap function
         gap0 = gap * he3_bcsgap_fast(0D0)/he3_bcsgap_fast(ttc)
!        k = (2.5D0 - gap0)
!     .    / (1D0 - dsqrt(2D0*const_pi*gap)
!     .             *dexp(-gap)*(1D0 + 3D0/8D0/gap))

        k = 2.388693D0

        he3_yosida0_fast =
     .     dsqrt(2D0*const_pi*gap/ttc) * dexp(-gap/ttc)
     .     * (1D0 + 3D0/8D0 * ttc/gap) * (gap-gap/ttc)
     .     * (1D0 - ttc**k)
     .   + dexp(gap - gap/ttc) * ttc**(k-0.5D0)
      end

! Eizel-1991 f.90
      function he3_yosida_par(ttc, gap)
        implicit none
        real*8 ttc,gap
        include 'he3.fh'
        he3_yosida_par = (
     .       2D0 * he3_yosida(ttc, gap,0D0)
     .     + 3D0 * he3_yosida(ttc, gap,2D0)
     .    )/5D0
      end

! Eizel-1991 f.90
      function he3_yosida_perp(ttc, gap)
        implicit none
        real*8 ttc,gap
        include 'he3.fh'
        he3_yosida_par = (
     .       4D0 * he3_yosida(ttc, gap,0D0)
     .     + 1D0 * he3_yosida(ttc, gap,2D0)
     .    )/5D0
      end

! B-phase Normal component density \rho_n^B/\rho_0
! VW book f.3.92
      function he3_rho_nb(ttc, p)
        implicit none
        real*8 ttc,p,f1s,gap,Y0
        include 'he3.fh'
        f1s = He3_f1s(p)
        gap = he3_bcsgap(ttc)
        Y0  = He3_yosida(ttc, gap, 0D0)
        he3_rho_nb = (3D0+f1s)*Y0/(3D0+f1s*Y0)
      end

! He3-B suseptibility chi_b / chi_0
! see VW book ch.10 p.449; ch2 p.90
! see Wheatley-75 f 3.7
! There is also additional term to 3*chi0
!  + 2/5 F2a (1-Y0)^2
      function He3_chi_b(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,f0a,gap,Y0
        f0a = He3_f0a(p)
        gap = he3_bcsgap(ttc)
        Y0  = He3_yosida(ttc, gap, 0D0)
        He3_chi_b =
     .    (1D0 + f0a) * (2D0+Y0) /
     .    (3D0 + f0a * (2D0+Y0))
      end

! B-phase Leggett freq
      function he3_nu_b(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap
        gap  = he3_trivgap(ttc,p) * const_kb * he3_tc(p)/1D3 ! mk->K
        he3_nu_b = dsqrt(3D0 / 8D0 / const_pi /
     .                   he3_chi_b(ttc,p)/he3_chi_n(p))
     .    * he3_gyro**2 * const_hbar * he3_2n0(p) / 4D0
     .    * gap * dlog(he3_tfeff(p)*const_kB/gap)
     .    * dsqrt(0.9574D0 + 0.3682D0*dexp(-p/6.9234D0)) ! fit to experimental data
      end


! Suseptibility [sgs] vs P [bar], T [mK] -- Old
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



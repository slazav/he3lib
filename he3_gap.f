! BCS gap / (kB Tc) for pure 3He-B, t = T / Tc
! Newton iteration based on a note by EVT & RH
! From ROTA texture library
! see: http://ltl.tkk.fi/research/theory/qc/bcsgap.html
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
        ynew = 1.7638D0*DSQRT(1D0-ttc)/const_2pi
        do while (DABS(dy) > 1.0D-8)
          y = ynew
          root = DSQRT((dble(m)*ttc)**2+y**2)
          g = DLOG((dble(m)*ttc+root) / (2D0*dble(m)))
     .        - (1D0/dble(m)**2 - dble(m)*(ttc/root)**3)/24D0
          dg = y/(root*(dble(m)*ttc+root))
     .        - dble(m)*ttc**3*y/(8D0*root**5)
          do n=1,m
            root=DSQRT((ttc*(dble(n)-0.5D0))**2 + y**2)
            g = g + 1D0/(dble(n)-0.5D0) - ttc/root
            dg = dg + ttc*y/root**3
          end do
          dy = g/dg
          ynew = ynew-dy
        end do
        he3_bcsgap = const_2pi*ynew
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

! Heat capacity jump for He3-B, DeltaCb/Cn
! Greywall-1986, Fig.19
      function he3_dcbcn(p)
        implicit none
        include 'he3.fh'
        real*8 p
        he3_dcbcn = 41.9D0 / he3_vm(p) + 0.322D0
      end

! Heat capacity jump for He3-A, DeltaCa/Cn
! Greywall-1986, Fig.19
      function he3_dcacn(p)
        implicit none
        include 'he3.fh'
        real*8 p
        he3_dcacn = 94.2D0 / he3_vm(p) - 1.58D0
      end


! Trivial strong-coupling correction to the BCS energy gap
! see: http://ltl.tkk.fi/research/theory/qc/bcsgap.html
! Corrections are from Serene,Rainer-1983 (Phys. Rep. 101, 221), table 4
! Specific heat jump (dcpcn) is from Greywall-1985 paper (PRB 33, 7520), Fig.19
!      function he3_trivgap(ttc,p)
!        implicit none
!        include 'he3.fh'
!        integer it, ic
!        real*8 ttc, p, dcpcn, wt1, wt2, wc1, wc2, corr
!        real*8 c,x
!        dimension c(5), x(11,5)
!        c = (/ 1.43D0,1.6D0,1.8D0,2.D0,2.2D0 /)
!        dcpcn = he3_dcbcn(p)
!        x(11,1:5) = (/ 1.D0,1.056D0,1.115D0,1.171D0,1.221D0 /)
!        x(10,1:5) = (/ 1D0,1.048D0,1.097D0,1.141D0,1.18D0 /)
!        x(9,1:5)  = (/ 1D0,1.041D0,1.083D0,1.119D0,1.15D0 /)
!        x(8,1:5)  = (/ 1D0,1.036D0,1.072D0,1.102D0,1.128D0 /)
!        x(7,1:5)  = (/ 1D0,1.032D0,1.063D0,1.089D0,1.112D0 /)
!        x(6,1:5)  = (/ 1D0,1.028D0,1.056D0,1.079D0,1.099D0 /)
!        x(5,1:5)  = (/ 1D0,1.026D0,1.051D0,1.073D0,1.091D0 /)
!        x(4,1:5)  = (/ 1D0,1.024D0,1.049D0,1.069D0,1.086D0 /)
!        x(3,1:5)  = (/ 1D0,1.024D0,1.048D0,1.068D0,1.085D0 /)
!        x(2,1:5)  = (/ 1D0,1.024D0,1.048D0,1.068D0,1.085D0 /)
!        x(1,1:5)  = (/ 1D0,1.024D0,1.048D0,1.068D0,1.085D0 /)
!        it=INT(ttc*10D0 - 1D-5) + 1
!        wt1 = (0.1D0*dble(it)-ttc)/0.1D0
!        wt2 = 1D0 - wt1
!        ic = 1
!        do
!          if (dcpcn < c(ic+1) ) exit
!          ic = ic+1
!          if (ic == 4) exit
!        end do
!        wc1 = (c(ic+1)-dcpcn)/(c(ic+1)-c(ic))
!        wc2 = 1D0 - wc1
!        corr = wt1*(wc1*x(it,ic)+wc2*x(it,ic+1))
!        corr = corr + wt2*(wc1*x(it+1,ic)+wc2*x(it+1,ic+1))
!        he3_trivgap = he3_bcsgap(ttc)*corr
!      end

! Trivial strong-coupling correction to the BCS energy gap.
! Application of Einzel-2003 interpolation formula to Serene,Rainer-1983 
! corrections (Phys. Rep. 101, 221), table 4
      function he3_trivgap(ttc,p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 y,l1,l2,g0,ga,aa,x
        real*8 yb,l1b,l2b,g0b,aab,gab  ! BCS values

        ! see my notes somewhere
        x=1D0-ttc
        y=he3_dcbcn(p)

        yb = 12D0/7D0/const_z3
        l1b = 2D0/3D0*yb
        l2b = yb*(1D0-31D0/144D0*const_z5*yb**2)
        g0b = const_pi/dexp(const_euler)
        aab = -l2b/l1b + 2D0*l1b/3D0 * dexp(2D0*const_euler) - 1D0
        gab = g0b*dtanh(const_pi/g0b*dsqrt(l1b/(1D0-x)*(x+aab*x**2)))

        l1 = 2D0/3D0*y*(1.021745D0+0.002692D0*y-0.012582D0*y**2)
        l2 = y*(-0.436318D0+0.733004D0*y-0.032478D0*y**2)
        g0 = g0b*(0.700558D0+0.275132*y-0.045649D0*y**2)
        aa = -l2/l1 + 2D0*l1/3D0 * dexp(2D0*const_euler) - 1D0
        ga = g0*dtanh(const_pi/g0*dsqrt(l1/(1D0-x)*(x+aa*x**2)))

        he3_trivgap = ga/gab*he3_bcsgap(ttc)
      end


! delta(0)/Tc for 3He versus pressure, bar
! linear interpolation in density between BCS value at zero bar
! and Todoschenko's value 1.99 at melting pressure
      function he3_todogap(ttc,p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p, gap0, gap1
        real*8 r,r0,r1
        gap0 = he3_bcsgap(0D0)
        gap1 = 1.99D0
        r0 = 1D0/he3_vm(0D0)
        r1 = 1D0/he3_vm(34.338D0)
        r  = 1D0/he3_vm(p)
        he3_todogap = (gap0 + (r-r0)*(gap1-gap0)/(r1-r0)) *
     .                 he3_trivgap(ttc,p)/he3_trivgap(0D0,p)
      end

! wrapper function which should be used everywhere in the lib
      function he3_gap(ttc,p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        !! select one of gap functions
        !he3_gap = he3_bscgap(ttc)
        !he3_gap = he3_bcsgap_fast(ttc)
        he3_gap = he3_trivgap(ttc,p)
        !he3_gap = he3_todogap(ttc,p)
      end

! he3_gap expressed in energy units [erg] rather then $T_c$
      function he3_egap(ttc,p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        he3_egap = he3_gap(ttc,p)*he3_tc(p)/1D3*const_kb
      end




! Integrand for Yosida function calculations
! x = tanh(\xi)/2 change is made to get good integrand
! and [0:1] integrating range.  d\xi -> 2 dx / (1-x**2)
! See also tests/plot_yosida_int.m
      function he3_yosida_int(x)
        implicit none
        include 'he3.fh'
        real*8 x, he3_yosida_int
        real*8 ttc, gap
        common /he3_yosida_int_cb/ ttc, gap, n
        real*8 xi, ek, n, C
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
        include 'he3_math.fh'
        real*8 he3_yosida_int
        external he3_yosida_int
        real*8 ttc, gap, n
        real*8 ttc1, gap1, n1
        common /he3_yosida_int_cb/ ttc1, gap1, n1
        ttc1=ttc
        gap1=gap
        n1=n

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

        he3_yosida = math_dint(he3_yosida_int, 0D0, 1D0, 100, 0)
      end

! Integrand for Entropy Yosida function (see Einzel-2003, table 1)
! x = tanh(\xi)/2 change is made to get good integrand
! and [0:1] integrating range.  d\xi -> 2 dx / (1-x**2)
! See also tests/plot_yosida_int.m
      function he3_yosida_s_int(x)
        implicit none
        include 'he3.fh'
        real*8 x, he3_yosida_s_int
        real*8 ttc, gap
        common /he3_yosida_s_int_cb/ ttc, gap
        real*8 xi, ek, C
        C=4D0
        xi = datanh(x)*C
        ek=dsqrt(xi**2 + gap**2)
        he3_yosida_s_int =
     .     (xi/ttc)**2
     .     / (dcosh(ek/(2D0*ttc)))**2 / 2D0/ttc
     .     * C / (1D0-x**2)
      end


! Entropy Yosida function vs T/Tc, gap
! See Einzel-2004
      function he3_yosida_s(ttc, gap)
        implicit none
        include 'he3.fh'
        include 'he3_math.fh'
        real*8 he3_yosida_s_int
        external he3_yosida_s_int
        real*8 ttc, gap
        real*8 ttc1, gap1
        common /he3_yosida_s_int_cb/ ttc1, gap1
        ttc1=ttc
        gap1=gap

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

        he3_yosida_s = 3D0/const_pi**2
     .              * math_dint(he3_yosida_s_int, 0D0, 1D0, 100, 0)
      end

! Heat Capacity Yosida function vs T/Tc, P
! See D.Einzel-2003
! Just simple derivative of he3_yosida_s
! Note: argument P instead of gap is needed
      function he3_yosida_c(ttc, P)
        implicit none
        include 'he3.fh'
        include 'he3_math.fh'
        real*8 ttc, P
        real*8 dtp,dtm, gap,gapm,gapp
        dtp=ttc*1D-3
        dtm=ttc*1D-3

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
        if (ttc.gt.1D0-dtp) then
          dtp = 1D0-ttc
        endif

        if (P.lt.0D0) then
          gap = he3_bcsgap(ttc);
          gapp = he3_bcsgap(ttc+dtp);
          gapm = he3_bcsgap(ttc-dtm);
        else
          gap = he3_gap(ttc, P);
          gapp = he3_gap(ttc+dtp, P);
          gapm = he3_gap(ttc-dtm, P);
        endif

        ! Yc = C/T = dS/dT = d/dT (T*Ys) = T dYs/dT + Ys
        he3_yosida_c = ttc *
     .   (he3_yosida_s(ttc+dtp, gapp)-he3_yosida_s(ttc-dtm, gapm))/
     .   (dtp+dtm)
     .   + he3_yosida_s(ttc, gap)
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
        he3_yosida_perp = (
     .       4D0 * he3_yosida(ttc, gap,0D0)
     .     + 1D0 * he3_yosida(ttc, gap,2D0)
     .    )/5D0
      end

! Interpolation of Yosida function Y0 (Einzel-2003)
! See my notes somewhere
! If P<0 BCS is used
      function he3_yosida_0_appr(ttc,p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p 
        real*8 gap, gap_bcs, dcbcn, dcbcn_bcs
        real*8 y0, y0tc, dy0tc, k, s, bet

        gap_bcs = const_pi/dexp(const_euler)
        dcbcn_bcs = 12D0/7D0/const_z3
        if (P.lt.0D0) then
          gap  = gap_bcs
          dcbcn = dcbcn_bcs
        else
          gap = he3_trivgap(0D0,P) ! zero-temperature gap
          dcbcn = he3_dcbcn(P)
        endif

        bet = 0D0 ! 3D0/8D0
        s   = 2D0*dcbcn/dcbcn_bcs
        y0  = dsqrt(const_2pi*gap/ttc)*dexp(-gap/ttc)
     .      * (1D0 + bet*ttc/gap)
        y0tc = dsqrt(const_2pi*gap)*dexp(-gap)*(1D0+bet/gap)
        dy0tc = y0tc*(gap-0.5D0*(gap-bet)/(gap+bet))

        k = (s - dy0tc/y0tc)/(1D0-y0tc);
        he3_yosida_0_appr = y0*(1D0-(1D0-1D0/y0tc)*ttc**k)
      end

! Interpolation of Entropy Yosida function Ys (Einzel-2003)
      function he3_yosida_s_appr(ttc,p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p 
        real*8 gap, gap_bcs, dcbcn, dcbcn_bcs
        real*8 y0, y0tc, dy0tc, k, s, bet

        gap_bcs = const_pi/dexp(const_euler)
        dcbcn_bcs = 12D0/7D0/const_z3
        if (P.lt.0D0) then
          gap  = gap_bcs
          dcbcn = dcbcn_bcs
        else
          gap = he3_trivgap(0D0,P) ! zero-temperature gap
          dcbcn = he3_dcbcn(P)
        endif

        bet = 0D0 ! 15D0/8D0
        s   = dcbcn
        y0   = 3D0/const_pi**2 * gap/ttc
     .       * dsqrt(const_2pi*gap/ttc)*dexp(-gap/ttc)
     .       * (1D0 + bet*ttc/gap)
        y0tc = 3D0/const_pi**2 * gap
     .       * dsqrt(const_2pi*gap)*dexp(-gap)
     .       * (1D0 + bet/gap)
        dy0tc = (gap - 1.5D0 + bet/(gap+bet))*y0tc

        k = (s - dy0tc/y0tc)/(1D0-y0tc);
        he3_yosida_s_appr = y0*(1D0-(1D0-1D0/y0tc)*ttc**k)
      end

! Interpolation of Heat capacity Yosida function Yc (Einzel-2003)
      function he3_yosida_c_appr(ttc,p)
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 gap, gap_bcs, dcbcn, dcbcn_bcs
        real*8 y0, y0tc, dy0tc, ytc, k, s, bet, l2

        gap_bcs = const_pi/dexp(const_euler)
        dcbcn_bcs = 12D0/7D0/const_z3
        if (P.lt.0D0) then
          gap  = gap_bcs
          dcbcn = dcbcn_bcs
          l2 = 0.77865D0
        else
          gap = he3_trivgap(0D0,P) ! zero-temperature gap
          dcbcn = he3_dcbcn(P)
          l2 = dcbcn *
     .      (-0.032478D0*dcbcn**2D0 + 0.733004D0*dcbcn - 0.436318D0)
        endif

        bet = 0D0 ! 11D0/8D0
        ytc = 1D0 + dcbcn
        ! extracted from S-R gap corrections
        s = 3D0*l2/ytc

        y0   = 3D0*(gap/ttc/const_pi)**2
     .       * dsqrt(const_2pi*gap/ttc)*dexp(-gap/ttc)
     .       * (1D0 + bet*ttc/gap)
        y0tc = 3D0*(gap/const_pi)**2
     .       * dsqrt(const_2pi*gap)*dexp(-gap)
     .       * (1D0 + bet/gap)
        dy0tc = (gap - 2.5D0 + bet/(gap+bet))*y0tc

        k = (s - dy0tc/y0tc)/(1D0-y0tc/ytc)
        he3_yosida_c_appr = y0*(1D0-(1D0-ytc/y0tc)*ttc**k)
      end

! Z3,Z5,Z7, lambda
! Code from http://ltl.tkk.fi/research/theory/qc/bcsgap.html
! Original nsplit=10 is too small for (z3 - 0.9*z5 + 0.9*z5.^2./z3 - 1.5*z7)
! combination in he3_text_lhv
      function he3_z3(ttc,gap)
        implicit none
        include 'he3.fh'
        real*8 ttc,gap
        real*8 x,xs,tm,sq
        integer i, nsplit
        parameter (nsplit=100)

        if (ttc.lt.0D0) then
          he3_z3=NaN
          return
        endif

        tm=ttc*dfloat(nsplit)
        x=gap/const_2pi
        xs=x**2
        sq=dsqrt(tm**2+xs)
        he3_z3 = 1D0/(sq*(tm+sq)) - tm*ttc**2/( 8D0*sq**5)
        do i=1,nsplit
          sq=dsqrt((ttc*(dfloat(i)-0.5D0))**2+xs)
          he3_z3 = he3_z3 + ttc/sq**3
        enddo
        he3_z3 = he3_z3 / const_2pi**2 * gap*gap
      end

      function he3_z5(ttc,gap)
        implicit none
        include 'he3.fh'
        real*8 ttc,gap
        real*8 x,xs,tm,sq
        integer i, nsplit
        parameter (nsplit=100)

        if (ttc.lt.0D0) then
          he3_z5=NaN
          return
        endif

        tm=ttc*dfloat(nsplit)
        x=gap/const_2pi
        xs=x**2
        sq=dsqrt(tm**2+xs)
        he3_z5 = (tm+2D0*sq)/(3D0*sq**3*(tm+sq)**2)
     .    - 5D0*tm*ttc**2/(24D0*sq**7)
        do i=1,nsplit
          sq=dsqrt((ttc*(dfloat(i)-0.5D0))**2+xs)
          he3_z5 = he3_z5 + ttc/sq**5
        enddo
        he3_z5 = he3_z5 * xs/const_2pi**2 * gap*gap
      end

      function he3_z7(ttc,gap)
        implicit none
        include 'he3.fh'
        real*8 ttc,gap
        real*8 x,xs,tm,sq
        integer i, nsplit
        parameter (nsplit=100)

        if (ttc.lt.0D0) then
          he3_z7=NaN
          return
        endif

        tm=ttc*dfloat(nsplit)
        x=gap/const_2pi
        xs=x**2
        sq=dsqrt(tm**2+xs)
        he3_z7 = (11D0*tm**2+9D0*tm*sq+8D0*xs)
     .            / (15D0*sq**5*(tm+sq)**3)
     .           - 7D0*tm*ttc**2/(24D0*sq**9)
        do i=1,nsplit
          sq=dsqrt((ttc*(dfloat(i)-0.5D0))**2+xs)
          he3_z7 = he3_z7 + ttc/sq**7
        enddo
        he3_z7 = he3_z7 * xs**2/const_2pi**2 * gap*gap
      end

      function he3_lambda(ttc,gap)
        implicit none
        include 'he3.fh'
        real*8 ttc,gap
        real*8 x,xs,tm,sq
        integer i, nsplit
        parameter (nsplit=100)

        if (ttc.lt.0D0) then
          he3_lambda=NaN
          return
        endif

        tm=ttc*dfloat(nsplit)
        x=gap/const_2pi
        xs=x**2
        sq=dsqrt(tm**2+xs)
        he3_lambda = 1D0 - tm/(x+sq) -
     .     tm*ttc**2*x*(x+2D0*sq)/(24D0*sq**3*(x+sq)**2)
        do i=1,nsplit
          sq=dsqrt((ttc*(dfloat(i)-0.5D0))**2+xs)
          he3_lambda = he3_lambda + ttc*x/(sq*(sq+x))
        enddo
      end

! B-phase Normal component density \rho_n^B/\rho_0
! VW book f.3.92
      function he3_rho_nb(ttc, p)
        implicit none
        real*8 ttc,p,gap,f1s,Y0
        include 'he3.fh'
        f1s = He3_f1s(p)
        gap = he3_gap(ttc,p)
        Y0  = He3_yosida(ttc, gap, 0D0)
        he3_rho_nb = (3D0+f1s)*Y0/(3D0+f1s*Y0)
      end

! He3-B susceptibility chi_b / chi_0
! see VW book ch.10 p.449; ch2 p.90
! see Wheatley-75 f 3.7
! There is also additional term to 3*chi0
!  + 2/5 F2a (1-Y0)^2
      function He3_chi_b(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap,f0a,Y0
        f0a = He3_f0a(p)
        gap = he3_gap(ttc,p)
        Y0  = He3_yosida(ttc, gap, 0D0)
        He3_chi_b =
     .    (1D0 + f0a) * (2D0+Y0) /
     .    (3D0 + f0a * (2D0+Y0))
      end

! He3-B Cooper pair susceptibility ratio chi_bp / chi_b
! see Leggett-Takagi 1975, f.12
      function He3_chi_bp(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap,Y0,Y2
        gap = he3_gap(ttc,p)
        Y0  = He3_yosida(ttc, gap, 0D0)
        Y2  = He3_yosida(ttc, gap, 2D0)
        He3_chi_bp =
     .    2D0*(1D0 - Y2) / (2D0+Y0)
      end

! He3-B heat capacity (C/R)
      function he3_c_b(ttc,P)
        implicit none
        include 'he3.fh'
        real*8 ttc,P
        he3_c_b = he3_gammaf(P)
     .          * ttc*he3_tc(P)*1D-3   ! T,K
     .          * he3_yosida_c(ttc,P)
      end


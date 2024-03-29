!HH> BCS gap + strong coupling corrections for He3-B

!> BCS energy gap, $\Delta/k_BT_c$ vs $T/T_c$
!> Newton iteration based on a note by E.Thuneberg and R.Hanninen.
!> <br>Taken from ROTA texture library.
!> See: <a href="http://ltl.tkk.fi/research/theory/qc/bcsgap.html">[1]</a>
!>      <a href="http://ltl.tkk.fi/research/theory/qc/bcsgap.pdf">[2]</a>
      function he3_bcsgap(ttc) !F>
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

!> BCS energy gap, $\Delta/k_BT_c$ vs $T/T_c$, Einzel approximation
!> D.Einzel JLTP 84 (1991) f.68.
!> <0.5% accuracy in the whole temperature range,
!> 70 times faster then he3_bcsgap
      function he3_bcsgap_fast(ttc) !F>
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

!> Derivative of BCS energy gap, $d\Delta^2/d(T/Tc)$ vs $T/T_c$
!> Same method as in he3_bcsgap calculation, V.Zavjalov, 2020
      function he3_bcsdgap2(ttc) !F>
        implicit none
        include 'he3.fh'
        integer n, m
        real*8 ttc,gap, root,y,dy,ynew,g1,g2
        gap = he3_bcsgap(ttc)/const_2pi
        if (ttc.eq.1D0) then
          ! See Einzel-2003, Eq 3
          he3_bcsdgap2 = -8D0/7D0/const_z3 * const_pi**2
          return
        end if
        if (ttc.gt.1D0) then
          he3_bcsdgap2 = 0D0
          return
        end if
        if (ttc.lt.0D0) then
          he3_bcsdgap2 = NaN
          return
        end if
        m = 30
        ! From the gap formula:
        !  dg/dt + dg/dy^2 * dy^2/dt = 0

        root = DSQRT((dble(m)*ttc)**2+gap**2)
        g1 = (dble(m) + dble(m)**2*ttc/root)/(dble(m)*ttc+root)
     .     + dble(m)*ttc**2/root**3/8D0
     .     - dble(m)*ttc**3/root**5/8D0 * dble(m)**2 * ttc

        g2 = 0.5D0 /root / (dble(m)*ttc+root)
     .     - 0.5D0 * dble(m)*ttc**3/root**5/8D0
        do n=1,m
          root=DSQRT((ttc*(dble(n)-0.5D0))**2 + gap**2)
          g1 = g1 - 1D0/root
     .           + ttc/root**3 * ttc*(dble(n)-0.5D0)**2
          g2 = g2 + 0.5D0 * ttc/root**3
        end do
        he3_bcsdgap2 = -g1/g2 * const_2pi**2
      end

!> Heat capacity jump for He3-B, $\Delta C_b/C_n$ vs P [bar], (exp.data, Greywall-1986)
!> Greywall-1986, Fig.19
      function he3_dcbcn(p) !F>
        implicit none
        include 'he3.fh'
        real*8 p
        he3_dcbcn = 41.9D0 / he3_vm(p) + 0.322D0
      end

!> Heat capacity jump for He3-A, $Delta C_a/C_n$ vs P [bar], (exp.data, Greywall-1986)
!> Greywall-1986, Fig.19
      function he3_dcacn(p) !F>
        implicit none
        include 'he3.fh'
        real*8 p
        he3_dcacn = 94.2D0 / he3_vm(p) - 1.58D0
      end

!> Trivial strong-coupling correction (WCP) to the BCS energy gap. $\Delta/k_BT_c$ vs $T/T_c$, P[bar]
!> Approximation of Serene,Rainer-1983 corrections (Phys.Rep. 101, 221), table 4.
!> Note that derivative of
!> the $\Delta^2$ in $T_c$ is not strictly proportional to the heat capacity jump.
!> This shows that exact heat capacity calculation requires WCP energy terms,
!> not just BCS calculations with modified gap...
      function he3_trivgap(ttc,p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 k,x,y,a,b,c,d

        y = he3_dcbcn(p) - 12D0/7D0/const_z3
        x = 1D0-ttc
        a = -0.08342397526067D0*y**2
     .     + 0.28724078111814D0*y + 1D0
        b =  0.03832526906988D0*y**2
     .     + 0.38097566888497D0*y
        c =  0.22188897390807D0*y**2
     .     + 1.33772124075848D0*y
     .     + 2.16590470748955D0
        d = -2.04760212642877D0*y**2
     .     + 0.59148230270250D0*y
     .     + 2.79334876756980D0

        k = a + b*dexp(-x*c-x**2*d)

        he3_trivgap = dsqrt(k)*he3_bcsgap(ttc)
      end

!> Derivative of the trivial strong-coupling (WCP) gap: d(Delta^2)/d(T/Tc)
      function he3_trivdgap2(ttc,p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        real*8 k,dk,x,y,a,b,c,d

        y = he3_dcbcn(p) - 12D0/7D0/const_z3
        x = 1D0-ttc
        a = -0.08342397526067D0*y**2
     .     + 0.28724078111814D0*y + 1D0
        b =  0.03832526906988D0*y**2
     .     + 0.38097566888497D0*y
        c =  0.22188897390807D0*y**2
     .     + 1.33772124075848D0*y
     .     + 2.16590470748955D0
        d = -2.04760212642877D0*y**2
     .     + 0.59148230270250D0*y
     .     + 2.79334876756980D0
        k = a + b*dexp(-x*c-x**2*d)
        dk = b*dexp(-x*c-x**2*d)*(c+2D0*x*d) ! d/dx = -d/dttc

        ! g^2 = g0^2 * k
        ! d(g^2) = d(g0^2)*k + g0^2 * dk

        he3_trivdgap2 = he3_bcsdgap2(ttc) * k
     .                + he3_bcsgap(ttc)**2 * dk
      end

!> Gap corrected to Todoschenko's value 1.99 at T=0,P=Pmelt, delta/Tc vs T/Tc, P[bar]
!> Linear interpolation in density between BCS value at zero bar
!> and Todoschenko's value 1.99 at melting pressure
      function he3_todogap(ttc,p) !F>
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


!> Wrapper function which should be used everywhere in the lib, same as he3_trivgap
      function he3_gap(ttc,p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        !! select one of gap functions
        !he3_gap = he3_bscgap(ttc)
        !he3_gap = he3_bcsgap_fast(ttc)
        he3_gap = he3_trivgap(ttc,p)
        !he3_gap = he3_todogap(ttc,p)
      end

!> he3_gap expressed in energy units [erg] rather then $T_c$
      function he3_egap(ttc,p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc, p
        he3_egap = he3_gap(ttc,p)*he3_tc(p)/1D3*const_kb
      end




! Integrand for Yosida function calculations
! x = tanh(\xi)/2 change is made to get good integrand
! and [0:1] integrating range.  d\xi -> 2 dx / (1-x**2)
! See also tests/he3_yosida
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


!> Yosida function of order n vs T/Tc, gap
!> See D.Einzel JLTP 84
!> <br>$Y_n = 2\int_0^\infty \left(\frac{\xi}{E}\right)^n\left(-\frac{\partial f^0}{\partial E}\right)$
!> <br>At T -> 0:
!> <br>$Y_n = 2\Gamma\left(\frac{n+1}{2}\right)\left(\frac{2k_BT}{\Delta}\right)^{(n-1)/2} \exp\left(-\frac{\Delta}{k_BT}\right)$
      function he3_yosida(ttc, gap, n) !F>
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
! See also tests/he3_yosida
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

!> Entropy Yosida function vs T/Tc, gap, see Einzel-2004
      function he3_yosida_s(ttc, gap) !F>
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
          he3_yosida_s=NaN
          return
        endif
        if (ttc.eq.0D0) then
          he3_yosida_s=0D0
          return
        endif
        if (ttc.gt.1D0) then
          he3_yosida_s=1D0
          return
        endif

        he3_yosida_s = 3D0/const_pi**2
     .              * math_dint(he3_yosida_s_int, 0D0, 1D0, 100, 0)
      end

! Integrand for Heat capacity Yosida function (see Einzel-2003, table 1)
! x = tanh(\xi)/2 change is made to get good integrand
! and [0:1] integrating range.  d\xi -> 2 dx / (1-x**2)
! See also tests/he3_yosida
      function he3_yosida_c_int(x)
        implicit none
        include 'he3.fh'
        real*8 x, he3_yosida_c_int
        real*8 ttc, gap, dgap2
        common /he3_yosida_c_int_cb/ ttc, gap, dgap2
        real*8 xi, ek, C
        C=4D0
        xi = datanh(x)*C
        ek=dsqrt(xi**2 + gap**2)
        he3_yosida_c_int =
     .     (Ek**2/ttc - 0.5D0*dgap2)/ttc
     .     / (dcosh(ek/(2D0*ttc)))**2 / 2D0/ttc
     .     * C / (1D0-x**2)
      end

!> Heat Capacity Yosida function vs T/Tc, gap, dgap2, see D.Einzel-2003
      function he3_yosida_c(ttc, gap, dgap2) !F>
        implicit none
        include 'he3.fh'
        include 'he3_math.fh'
        real*8 he3_yosida_c_int
        external he3_yosida_c_int
        real*8 ttc, gap, dgap2
        real*8 ttc1, gap1, dgap21
        common /he3_yosida_c_int_cb/ ttc1, gap1, dgap21
        ttc1=ttc
        gap1=gap
        dgap21=dgap2

        if (ttc.lt.0D0) then
          he3_yosida_c=NaN
          return
        endif
        if (ttc.eq.0D0) then
          he3_yosida_c=0D0
          return
        endif
        if (ttc.gt.1D0) then
          he3_yosida_c=1D0
          return
        endif

        he3_yosida_c = 3D0/const_pi**2
     .              * math_dint(he3_yosida_c_int, 0D0, 1D0, 100, 0)
      end


!> $Y_\parallel = 2/5 Y_0 + 3/5 Y_2$, see Eizel-1991 f.90
      function he3_yosida_par(ttc, gap) !F>
        implicit none
        real*8 ttc,gap
        include 'he3.fh'
        he3_yosida_par = (
     .       2D0 * he3_yosida(ttc, gap,0D0)
     .     + 3D0 * he3_yosida(ttc, gap,2D0)
     .    )/5D0
      end

!> $Y_\perp = 4/5 Y_0 + 1/5 Y_2$, see Eizel-1991 f.90
      function he3_yosida_perp(ttc, gap) !F>
        implicit none
        real*8 ttc,gap
        include 'he3.fh'
        he3_yosida_perp = (
     .       4D0 * he3_yosida(ttc, gap,0D0)
     .     + 1D0 * he3_yosida(ttc, gap,2D0)
     .    )/5D0
      end

!H> Z3,Z5,Z7, and lambda functions
!> Code from http://ltl.tkk.fi/research/theory/qc/bcsgap.html
!> Original nsplit=10 is too small for (z3 - 0.9*z5 + 0.9*z5.^2./z3 - 1.5*z7)
!> combination in he3_text_lhv

!> Z3 function
      function he3_z3(ttc,gap) !F>
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

!> Z5 function
      function he3_z5(ttc,gap) !F>
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

!> Z7 function
      function he3_z7(ttc,gap) !F>
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

!> Lambda function
      function he3_lambda(ttc,gap) !F>
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

!H> B-phase normal component density, susceptibility, and heat capacity

!> B-phase Normal component density \rho_n^B/\rho_0
!> VW book f.3.92
      function he3_rho_nb(ttc, p) !F>
        implicit none
        real*8 ttc,p,gap,f1s,Y0
        include 'he3.fh'
        f1s = he3_f1s(p)
        gap = he3_gap(ttc,p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        he3_rho_nb = (3D0+f1s)*Y0/(3D0+f1s*Y0)
      end

!> He3-B susceptibility chi_b/chi_0
!> see VW book ch.10 p.449, ch2 p.90;
!> see Wheatley-75 f 3.7;
!> There is also additional term to 3*chi0: + 2/5 F2a (1-Y0)^2
      function he3_chi_b(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap,f0a,Y0
        f0a = he3_f0a(p)
        gap = he3_gap(ttc,p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        he3_chi_b =
     .    (1D0 + f0a) * (2D0+Y0) /
     .    (3D0 + f0a * (2D0+Y0))
      end

!> He3-B Cooper pair susceptibility ratio chi_bp/chi_b
!> see Leggett-Takagi 1975, f.12
      function he3_chi_bp(ttc, p) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap,Y0,Y2
        gap = he3_gap(ttc,p)
        Y0  = he3_yosida(ttc, gap, 0D0)
        Y2  = he3_yosida(ttc, gap, 2D0)
        he3_chi_bp =
     .    2D0*(1D0 - Y2) / (2D0+Y0)
      end

!> He3-B heat capacity (C/R)
      function he3_c_b(ttc,P) !F>
        implicit none
        include 'he3.fh'
        real*8 ttc,P
        real*8 gap, dgap2
        gap = he3_trivgap(ttc, P)
        dgap2 = he3_trivdgap2(ttc, P)
        he3_c_b = he3_gammaf(P)
     .          * ttc*he3_tc(P)*1D-3   ! T,K
     .          * he3_yosida_c(ttc,gap,dgap2)
      end

!> <p>Some gap-related functions on the plot:
!> <p><img src="img/he3_gap.png">

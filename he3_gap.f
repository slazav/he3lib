! Yosida function vs T/Tc
! Sourse: B-phase notebook.
! Least squares fitting. From file: YOSHID
! Polinom of the order : 5
! Residual: 0.000
! Origin: Mukharskii, Dmitriev

      function He3_yosida(TTC)
        implicit none
        include 'he3.fh'
        real*8 TTC, A(5)
        real*8 XMIN,XMAX,XCAP
        integer IFAIL,M1
        DATA M1/5/
        DATA XMIN/9.9999994D-02/,XMAX/1.000000D0/
        DATA A/.7349111D0, .5123515D0, .1371038D0,
     .         -1.4855450D-02, -4.5979050D-03/
        if (TTC.GE.0.1D0) then
          IFAIL=1
          XCAP=((TTC-XMIN)-(XMAX-TTC))/(XMAX-XMIN)
          call E02AEE(M1,A,XCAP,He3_yosida,IFAIL)
          if (IFAIL.NE.0)print *,'Error in E02AEE :',IFAIL
        else
          He3_yosida = NaN
        end if
      end

! BCS gap / (kB Tc) for pure 3He-B, t = T / Tc
! Newton iteration based on a note by EVT & RH
! From ROTA texture library
      function he3_bcsgap(ttc)
        implicit none
        include 'he3.fh'
        integer n, m
        real*8 ttc,root,y,dy,ynew,g,dg
        m = 30
        dy = 1.0
        ynew = 1.7638*SQRT(1-ttc)/(2*const_pi)
        do while (ABS(dy) > 1.0E-8)
          y = ynew
          root = SQRT((m*ttc)**2+y**2)
          g = LOG((m*ttc+root)/(2*m))
     .        - (1D0/m**2-m*(ttc/root)**3)/24
          dg = y/(root*(m*ttc+root))
     .        - m*ttc**3*y/(8D0*root**5)
          DO n=1,m
            root=SQRT((ttc*(n-0.5D0))**2+y**2)
            g = g + 1D0/(n-0.5D0) - ttc/root
            dg = dg + ttc*y/root**3
          end do
          dy = g/dg
          ynew = ynew-dy
        end do
        he3_bcsgap = 2*const_pi*ynew
        if (ttc >= 1.0) he3_bcsgap = 0D0
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
        wt1 = (0.1D0*it-ttc)/0.1D0
        wt2 = 1D0 - wt1
        ic = 1
        do
          if (dcpcn < c(ic+1) ) exit
          ic = ic+1
          if (ic == 4) exit
        end do
        wc1 = (c(ic+1)-dcpcn)/(c(ic+1)-c(ic))
        wc2 = 1 - wc1
        corr = wt1*(wc1*x(it,ic)+wc2*x(it,ic+1))
        corr = corr + wt2*(wc1*x(it+1,ic)+wc2*x(it+1,ic+1))
        he3_trivgap = he3_bcsgap(ttc)*corr
      end

      function he3_z3(ttc,gap)
        implicit none
        include 'he3.fh'
        integer i,maxi
        real*8 ttc,gap,y,help,mt,corr1,corr2,sum
        y = gap/(2*const_pi)
        sum = 0D0
        maxi = 100
        mt = maxi*ttc
        do i=1,maxi
           sum = sum + ttc/((ttc*(i-0.5D0))**2+y**2)**1.5D0
        end do
        help = SQRT(mt**2+y**2)
        corr1 = 1/(help*(mt+help))
        corr2 = mt**3/(8*help**5)
        he3_z3 = y**2*(sum + corr1 - corr2)
      end

      function he3_z5(ttc,gap)
        implicit none
        include 'he3.fh'
        integer i,maxi
        real*8 ttc,gap,y,help,mt,corr1,corr2,sum
        y = gap/(2*const_pi)
        sum = 0D0
        maxi = 100
        mt = maxi*ttc
        do i=1,maxi
           sum = sum + ttc/((ttc*(i-0.5D0))**2 + y**2)**2.5D0
        end do
        help = SQRT(mt**2+y**2)
        corr1 = (mt+2*help)/(3*help**3*(mt+help)**2)
        corr2 = 5*mt**3/(24*help**7)
        he3_z5 = y**4*(sum + corr1 - corr2)
      end

      function he3_z7(ttc,gap)
        implicit none
        include 'he3.fh'
        integer i,maxi
        real*8 ttc,gap,y,help,mt,corr1,corr2,sum
        y = gap/(2*const_pi)
        sum = 0D0
        maxi = 100
        mt=maxi*ttc
        do i=1,maxi
           sum = sum + ttc/((ttc*(i-0.5D0))**2+y**2)**3.5D0
        end do
        help = SQRT(mt**2+y**2)
        corr1 = (11*mt*mt+9*mt*help+8*y*y)/(15*help**5*(mt+help)**3)
        corr2 = 7*mt**3/(24*help**9)
        he3_z7 = y**6*(sum + corr1 - corr2)
      end

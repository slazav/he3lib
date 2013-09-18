!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     heat capacity, Cv
!     Greywall-1983
      function He3_cv_n(t, v)
        implicit none
        include 'he3.fh'
        real*8 t,v,a,b,c,d
        real*8 s1,s2,s3
        integer i,j
        dimension a(5,4), b(4,3), c(3,3), d(3)
        a(1,1) = -2.9190414D0
        a(1,2) =  5.2893401D+2
        a(1,3) = -1.8869641D+4
        a(1,4) =  2.6031315D+5
        a(2,1) =  0D0
        a(2,2) =  0D0
        a(2,3) =  0D0
        a(2,4) =  0D0
        a(3,1) = -2.4752597D+3
        a(3,2) =  1.8377260D+5
        a(3,3) = -3.4946553D+6
        a(3,4) =  0D0
        a(4,1) =  3.8887481D+4
        a(4,2) = -2.8649769D+6
        a(4,3) =  5.2526785D+7
        a(4,4) =  0D0
        a(5,1) = -1.7505655D+5
        a(5,2) =  1.2809001D+7
        a(5,3) = -2.3037701D+8
        a(5,4) =  0D0

        b(1,1) = -6.5521193D-2
        b(1,2) =  1.3502371D-2
        b(1,3) =  0D0
        b(2,1) =  4.1359033D-2
        b(2,2) =  3.8233755D-4
        b(2,3) = -5.3468396D-5
        b(3,1) =  5.7976786D-3
        b(3,2) = -6.5611532D-4
        b(3,3) =  1.2689707D-5
        b(4,1) = -3.8374623D-4
        b(4,2) =  3.2072581D-5
        b(4,3) = -5.3038906D-7

        c(1,1) = -2.5482958D+1
        c(1,2) =  1.6416936D+0
        c(1,3) = -1.5110378D-2
        c(2,1) =  3.7882751D+1
        c(2,2) = -2.8769188D+0
        c(2,3) =  3.5751181D-2
        c(3,1) =  2.4412956D+1
        c(3,2) = -2.4244083D+0
        c(3,3) =  6.7775905D-2

        d(1) = -7.1613436D+0
        d(2) =  6.0525139D-1
        d(3) = -7.1295855D-3

        if (t.lt.0.1) then
          s1=0D0
          do i=1,5
            do j=0,3
              s1 = s1 + a(i,j+1) * v**(-j) * t**i
            enddo
          enddo
          he3_cv_n = s1
          return
        endif

        if (t.ge.0.1.and.t.lt.2.5) then
          s1=0D0
          do i=0,3
            do j=0,2
              s1 = s1 + b(i+1,j+1) * v**j * t**(-i)
            enddo
          enddo
          s2=0D0
          do i=1,3
            do j=0,2
              s2 = s2 + c(i,j+1) * v**j * t**(-i)
            enddo
          enddo
          s3=0D0
          do j=0,2
            s3 = s3 + d(j+1) * v**j
          enddo
          he3_cv_n = s1 + dexp(-s3/t) * s2
          return
        endif
        he3_cv_n = NaN
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!     .    * dsqrt(0.9574D0 + 0.3682D0*dexp(-p/6.9234D0)) ! fit to experimental data
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Osheroff's spin wave velocity. S [cm/s] vs P [bar], T [mK]
! Osheroff's spin wave velocity if recalculated to arbitrary pressure
! by taking value of Osheroff for melting pressure (1100 cm/sek)
! and assuming S is proportional to Fermi velocity.
! Origin: Mukharskii, Dmitriev

      function He3_swvel(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T
        He3_swvel =
     .   1100D0/He3_Vf(34.3D0)*He3_Vf(P)*SQRT(1D0-T/He3_Tc(P))
      end

! Parallel Fomin spin wave velocity. Cpar [cm/c] vs P [bar], T [mK]
! Origin: Mukharskii, Dmitriev

      function he3_swvel_par(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T
        he3_swvel_par=He3_swvel(P,T)*SQRT(2.0D0)
      end

! Perp Fomin spin wave velocity. Cper [cm/c] vs P [bar], T [mK]
! Origin: Mukharskii, Dmitriev

      function he3_swvel_per(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T
        he3_swvel_per=He3_swvel(P,T)*SQRT(1.5D0)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Leggett-Takagi tau_r and tau_s [s] vs T/Tc, 20bar
! Ref: WV pic.10.5 20bar

      function He3_tau_r(TTC)
        implicit none
        include 'he3.fh'
        real*8 TTC
        He3_tau_r = 1.2D-7/dsqrt(1.0D0-TTC)
      end

      function He3_tau_f(TTC)
        implicit none
        include 'he3.fh'
        real*8 TTC, P
        P=20D0
        He3_tau_f = 1D0 /
     .    (4D0*const_pi**2 *He3_nu_b(P,TTC) * He3_tau_r(TTC))
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!SPIN DIFFUSION COEFF, NORMAL D [cm**2/sec] vs P [bar], T[mk]
! Brewer data for D*T**2
! JLTP v.56 #5/6 p.617
! Least squares fitting. From file: YDIFD
! DT2=                  vs.  P=
! Polinom of the order : 4
! Residual: 0.000
! Origin: Mukharskii, Dmitriev

      function He3_Dn_exp(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T,DT2F, XCAP,F
        real*8 XMAX, XMIN
        real*8 A(4)
        integer IFAIL,M1
        DATA M1/4/
        DATA XMIN/0.0D0/,XMAX/27.75999D0/
        DATA A/1.0603673D-06, -5.4074212D-07,
     .         2.2418920D-07, -6.4375400D-08/
        IFAIL=1
        XCAP=((P-XMIN)-(XMAX-P))/(XMAX-XMIN)
        call E02AEE(M1,A,XCAP,F,IFAIL)
        if (IFAIL.NE.0) then
          print *,'Error in E02AEE :',IFAIL
          He3_Dn_exp = NaN;
        else if (he3_pmelt(T/1D3).lt.P.or.T.lt.He3_Tc(P)) then
          He3_Dn_exp = NaN;
        else
C       .89 accounts fo Grewall scale.
          DT2F=(F*0.89D0**2)
          He3_Dn_exp=DT2F/(T*1D-3)**2
        endif
      end

! SPIN DIFFUSION COEFF, SUPERFLUID D [cm**2/sec] vs P [bar], T[mk]
! Our data. Was measured at pressures .6,11,14.8,20,29 bar
! at T .45-.67,  .45-.8  .45-.81  .45-.94  .45-.67 Tc respectivly.
! Least squares 2-D fitting. From fil:: MUHLR1
! D/C= vs. T= & P=
! Polinoms of the orders : 4,  4
! Origin: Mukharskii, Dmitriev

      function He3_Ds_exp(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T,TR(1)
        real*8 WORK(5),A(25)
        real*8 XMIN,XMAX,YMIN,YMAX,F(1)
        integer IFAIL,K,L,NA,NWORK
        DATA XMIN/0.440D0/,XMAX/0.9400001D0/
        DATA YMIN/0.000D0/,YMAX/29.00000D0/
        DATA K/4/, L/4/, NA/25/, NWORK/5/
        DATA A/  6.6012963D-03, -2.2342042D-03,  4.2030680D-04
     ,  ,  1.6044091D-04, -8.2734930D-04,  1.7867290D-03, -6.6435340D-04
     ,  , -1.8627312D-04,  4.3953900D-04, -6.3450170D-04, -4.1724560D-04
     ,  ,  3.9544643D-04, -3.2521144D-04,  1.7599921D-04, -1.9650970D-04
     ,  ,  1.2244173D-04,  3.2786340D-05,  6.5047061D-08, -3.2710520D-05
     ,  , -6.1170401D-05,  3.7204154D-05,  9.9630890D-06,  1.2296370D-08
     ,  , -9.9427400D-06, -1.8589210D-05/
        TR(1)=T/He3_Tc(P)
        IFAIL=1
        call E02CBE(1,1,K,L,TR,XMIN,XMAX,P,YMIN,YMAX,F,
     ,              A,NA,WORK,NWORK,IFAIL)
        if (IFAIL.EQ.2) THEN
          print *,'Y out of range.'
          He3_Ds_exp=NaN
        else if (IFAIL.EQ.3) THEN
          print *,'X out of range.'
          He3_Ds_exp=NaN
        else if (IFAIL.NE.0) THEN
          print *,'Error:',IFAIL
          He3_Ds_exp=NaN
        else
          He3_Ds_exp = F(1)*he3_swvel_par(P,T)**(.6666666666666D0)
        end if
      end

! Spin diffusion coeff D [cm**2/sec] vs P [bar], T[mk]
! Origin: Mukharskii, Dmitriev

      function He3_D_exp(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T
        if (T.GT.He3_Tc(P)) THEN
          He3_D_exp = He3_Dn_exp(P,T)
        else
          He3_D_exp = He3_Ds_exp(P,T)
C          print *,'Superflow region. Check if data out range.'
        end if
      end

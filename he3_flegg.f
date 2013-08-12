! B-phase Leggett freq -- corrected to fit experimental data
! See /rota/Analysis/T-calibration/NMR/201308
! Correction is close to 2+F0a?
      function he3_exp_nu_b(ttc, p)
        implicit none
        include 'he3.fh'
        real*8 ttc,p,gap, corr
        corr = dsqrt(0.952D0 +0.6383D0*dexp(-p/6.627D0)) ! correction
        gap  = he3_trivgap(ttc,p) * const_kb * he3_tc(p)/1D3 ! mk->K
        ! same as in he3_nu_b
        he3_exp_nu_b = dsqrt(3D0 / 8D0 / const_pi / he3_chi_b(ttc,p))
     .    * he3_gyro**2 * const_hbar * he3_2n0(p) / 4D0
     .    * gap * dlog(he3_tfeff(p)*const_kB/gap) * corr
      end

! Legget freq^2, Lf^2 [Hz^2] vs P [bar], T/Tc
! Ahonen (18.7,21.1,25.4,29.0,32 bars), Osheroff(MP).
! Ahonen et.al.  JLTP. v.25 p.421(1976)
! Osheroff       PRL v.34. p.190
! From 0 to 18.7 bar data are not reliable.
! See file YOM3.
! Least squares 2-D fitting. From fil:: YOM3
! F= vs. T= & P=
! Polinmms of the orders : 4,  1
! Origin: Mukharskii, Dmitriev

      function He3_Flegg(P,T)
        implicit none
        include 'he3.fh'
        real*8 P,T,TT(1)
        real*8 WORK(5),A(10)
        real*8 XMIN,XMAX,YMIN,YMAX,F(1)
        integer IFAIL, K, L, NA, NWORK
        DATA XMIN/0.266D0/,XMAX/1.000D0/
        DATA YMIN/0.000D0/,YMAX/34.36D0/
        DATA K/4/,L/1/,NA/10/,NWORK/5/
        DATA A/  1.3169660D+11,  5.0310062D+10,
     .          -6.6371420D+10, -2.0536020D+10,
     .          -5.2457041D+09, -5.1179683D+09,
     .           5.8596372D+09,  3.1320310D+08,
     .          -6.9788020D+07,  2.0270360D+08/
        if (T.LT.He3_Tab(P)/He3_Tc(P)) THEN
          IFAIL=1
          TT(1)=T
          call E02CBE(1,1,K,L,TT,XMIN,XMAX,P,YMIN,YMAX,F,
     ,                A,NA,WORK,NWORK,IFAIL)
          He3_Flegg=F(1)
          if (IFAIL.EQ.2) THEN
!            print *,'P out of range:',YMIN,YMAX
            He3_Flegg=NaN
          elseif (IFAIL.EQ.3) THEN
!            print *,'T out of range:',XMIN,XMAX
            He3_Flegg=NaN
          elseif (IFAIL.NE.0) THEN
!            print *,'Error:',IFAIL
            He3_Flegg=NaN
          endif
          if (P.LT.18.7D0) print *,'Data not reliable.'
        else
!          print *,'No data for A-phase.'
          He3_Flegg=NaN
        endif
        return
      end

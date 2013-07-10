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
          He3_Dn_exp = -1D0;
        else if (he3_pmelt(T).lt.P.or.T.lt.He3_Tc(P)) then
          He3_Dn_exp = -1D0;
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
          He3_Ds_exp=-1D0
        else if (IFAIL.EQ.3) THEN
          print *,'X out of range.'
          He3_Ds_exp=-1D0
        else if (IFAIL.NE.0) THEN
          print *,'Error:',IFAIL
          He3_Ds_exp=-1D0
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

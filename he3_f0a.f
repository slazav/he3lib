! F0A vs P [bar]
! Sourse - Greywall, Phys.Rev.  v.27  5  (1983)
! Least squares fitting. From file: YF0A
! Polinom of the order : 7
! Residual: 0.000
! Origin: Mukharskii, Dmitriev

      function He3_F0a(P)
        implicit none
        include 'he3.fh'
        real*8  P, A(7)
        real*8 XMIN,XMAX,XCAP
        integer IFAIL,M1
        DATA M1/7/
        DATA XMIN/0.000000D0/,XMAX/34.36000D0/
        DATA A/-1.489332D0, -2.3159460D-02,  1.3571171D-02,
     .         -4.2908173D-03,  1.4413130D-03, -1.1601811D-03,
     .          9.9658221D-04/
        IFAIL=1
        XCAP=((P-XMIN)-(XMAX-P))/(XMAX-XMIN)
        call E02AEE(M1,A,XCAP,He3_F0a,IFAIL)
        if (IFAIL.NE.0) print *,'Error in E02AEE :',IFAIL
      end

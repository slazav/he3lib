! Z0 vs P [bar]
! Sourse B-phase notebook.
! Least squares fitting. From file: YZ0
! Polinom of the order : 5
! Residual: 0.000
! Origin: Mukharskii, Dmitriev

      function He3_z0(P)
        implicit none
        include 'he3.fh'
        real*8 P, A(5)
        real*8 XMIN,XMAX,XCAP
        integer IFAIL,M1
        DATA M1/ 5/
        DATA XMIN/0.000000D0/,XMAX/34.36000D0/
        DATA A/-5.762472D0, -0.1136529D0, 5.5511940D-02,
     .         -1.7914600D-02, 4.0055060D-03/
        IFAIL=1
        XCAP=((P-XMIN)-(XMAX-P))/(XMAX-XMIN)
        call E02AEE(M1,A,XCAP,He3_z0,IFAIL)
        if (IFAIL.NE.0) print *,'Error in E02AEE :',IFAIL
      end

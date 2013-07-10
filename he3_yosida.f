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
          He3_yosida = -1D0
        end if
      end
